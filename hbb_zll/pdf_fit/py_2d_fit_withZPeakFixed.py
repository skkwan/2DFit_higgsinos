# Run in ROOT 6.38 (do not do cmsenv)

import ROOT
from ROOT import RooFit as RF
import numpy as np
from array import array

##################################################
# Helper functions for spline (used for the signal MET)
##################################################
# From https://gitlab.cern.ch/cms-l1-ad/coffea-dask-axol1tl-studies/-/blob/master/prepareDatacards.py?ref_type=heads
def prune_knots_in_tiny_tail(knot_x, knot_y, tail_thresh=1e-8, min_keep=10):
    """
    Reduce knot density in the very low-yield tail without truncating the domain.
    Keeps first knot, last knot, and enough interior knots.

    tail_thresh: if remaining cumulative knot_y fraction is below this, stop adding dense knots.
    """
    knot_x = np.asarray(knot_x, dtype=np.float64)
    knot_y = np.asarray(knot_y, dtype=np.float64)

    if knot_x.size != knot_y.size or knot_x.size < 2:
        return knot_x, knot_y

    # Treat knot_y as "weights" proxy; compute cumulative from left
    w = np.maximum(knot_y, 0.0)
    total = float(np.sum(w))
    if total <= 0:
        return knot_x, knot_y

    keep = [0]
    cum = 0.0
    for i in range(1, len(knot_x) - 1):
        cum += w[i]
        remaining = total - cum
        # keep adding knots until remaining tail is tiny,
        # then stop keeping every knot (we'll keep last one).
        keep.append(i)
        if remaining / total < tail_thresh and len(keep) >= min_keep:
            break

    keep.append(len(knot_x) - 1)
    keep = np.unique(keep)
    return knot_x[keep], knot_y[keep]


def make_knot_x(
        x_min,
        x_max,
        n_knots,
        power = 1.0,
        min_dx_bins=2,        # enforce some minimum spacing in bin units
        centers=None,         # pass bin centers for snapping / spacing checks
    ):
        n_knots = int(max(2, n_knots))
        x_min = float(x_min)
        x_max = float(x_max)

        # Guard against x_min <= 0
        xmin_eff = max(x_min, 1e-3)
        #knot_x = np.exp(np.linspace(np.log(xmin_eff), np.log(x_max), n_knots))

        t = np.linspace(0.0, 1.0, n_knots)
        knot_x = x_min + (x_max - x_min) * (t**power)

        # Optional: snap to nearest bin center to avoid "between-bin" weirdness
        if centers is not None and len(centers) > 0:
            knot_x = np.array([centers[np.argmin(np.abs(centers - xx))] for xx in knot_x], dtype=float)

        # Optional: enforce minimum spacing (in bins) to avoid duplicates after snapping
        if centers is not None and min_dx_bins is not None and min_dx_bins > 0:
            # approximate bin width from centers
            if len(centers) > 1:
                bw = float(np.median(np.diff(centers)))
            else:
                bw = 1.0
            min_dx = min_dx_bins * bw

            filtered = [knot_x[0]]
            for xx in knot_x[1:]:
                if xx - filtered[-1] >= min_dx:
                    filtered.append(xx)
            if filtered[-1] != knot_x[-1]:
                filtered.append(knot_x[-1])
            knot_x = np.array(filtered, dtype=float)

        return knot_x

##################################################
##### Define fit observables 
##################################################
mll = ROOT.RooRealVar("m_ll", "m_ll", 60, 120)
met = ROOT.RooRealVar("met", "met", 0, 1200)

##################################################
###### Retrieve signal dataset (prepared with reformat.py)
##################################################
sigfilepath = 'snapshot_TChiZH_650_1_SR_mll_MET_fit_scheme.root'
sigfile = ROOT.TFile.Open(sigfilepath, "READ")
sigtree = sigfile.Get("event_tree")
weightXyear = ROOT.RooRealVar("weight_nominal", "weight_nominal", -1, 1)
variables = ROOT.RooArgSet(mll, met, weightXyear)
sigdataset = ROOT.RooDataSet("sigdataset", "sigdataset", variables, ROOT.RooFit.Import(sigtree), ROOT.RooFit.WeightVar(weightXyear))

##################################################
###### Retrieve background dataset (prepared with reformat.py)
##################################################
bkgfilepath = 'backgrounds_for_2D_fit.root'
bkgfile = ROOT.TFile.Open(bkgfilepath, "READ")
bkgtree = bkgfile.Get("event_tree")
weightXyear = ROOT.RooRealVar("weight_nominal", "weight_nominal", -1, 1)
variables = ROOT.RooArgSet(mll, met, weightXyear)
bkgdataset = ROOT.RooDataSet("bkgdataset", "bkgdataset", variables, ROOT.RooFit.Import(bkgtree), ROOT.RooFit.WeightVar(weightXyear))

###########################################################################
# Declare the signal 1D PDFs and multiply them into 2D PDFs
###########################################################################
# Get TH1 of signal met: https://github.com/guitargeek/hasco-2023-root/blob/main/notebooks/roofit-tutorial-01.ipynb 
sig_met_hist = ROOT.TH1D("sig_met_hist", "sig_met_hist", 120, 0, 1200)
sigtree.Draw("met >> sig_met_hist")
# Convert TH1 into RooDataHist
sig_roo_template_hist = ROOT.RooDataHist("sig_roo_template_hist", "sig_roo_template_hist", met, sig_met_hist)
# Create a RooHistPdf based on the RooFit histogram
sig_roohistpdf_met = ROOT.RooHistPdf("sig_roohistpdf_met", "sig_roohistpdf_met", met, sig_roo_template_hist, intOrder=0)

# Build spline-based signal met PDF from signal met histogram
nBins = sig_met_hist.GetNbinsX()
centers = np.array([sig_met_hist.GetBinCenter(i) for i in range(1, nBins+1)], dtype=np.float64)
vals    = np.array([max(float(sig_met_hist.GetBinContent(i)), 0.0) for i in range(1, nBins+1)], dtype=np.float64)

knot_avg_halfwidth_bins = 4
min_y = 1e-12
l=1.0
tail_thresh=1e-3
knot_x = make_knot_x(
    x_min=centers[0],
    x_max=centers[-1],
    n_knots=60,
    power = 2, #higher number = more dense at low mass
    centers=centers,
)
#find y value for each knot (average over local bins)
knot_y = []
for xx in knot_x:
    ib = int(np.argmin(np.abs(centers - xx)))
    i0 = max(0, ib - int(knot_avg_halfwidth_bins))
    i1 = min(nBins - 1, ib + int(knot_avg_halfwidth_bins))
    local = vals[i0:i1+1]
    yk = float(np.mean(local)) if local.size else float(vals[ib])
    knot_y.append(max(yk, min_y))


knot_x = np.array(knot_x, dtype=np.double)
knot_y = np.array(knot_y, dtype=np.double)

# From Claude: Thin sentinel floor values: keep a floor-valued knot only if it is directly
# adjacent to a non-floor knot. Removes long flat-zero plateaus that cause
# the cubic spline to oscillate negative 
is_floor = (knot_y <= min_y)
keep = np.array([
    i for i in range(len(knot_x))
    if not is_floor[i]
    or (i > 0 and not is_floor[i - 1])  # check point to the left (if we are not leftmost point), is it <floor?
    or (i < len(knot_x) - 1 and not is_floor[i + 1]) # repeat for the point to the right
])
knot_x = knot_x[keep]
knot_y = knot_y[keep]

# knot_x, knot_y = prune_knots_in_tiny_tail(knot_x, knot_y, tail_thresh=tail_thresh)
# knot_y_use = np.maximum(knot_y, min_y).astype(np.double)
vx = ROOT.std.vector('double')(knot_x)
vy_use = ROOT.std.vector('double')(knot_y)

print(list(zip(vx, vy_use)))

spline = ROOT.RooSpline(
        "spline", "spline",
        met,
        vx,
        vy_use,
        order=3
    )

pdf_of_spline = ROOT.RooGenericPdf(
    "pdf_of_spline", "pdf_of_spline",
    "max(@0, 1e-12)",    # prevent the interpolation from going 0 or negative
    ROOT.RooArgList(spline)
)

# # Also keep a RooHistPdf for reference
# sig_roo_template_hist = ROOT.RooDataHist("sig_roo_template_hist", "sig_roo_template_hist", ROOT.RooArgSet(met), sig_met_hist)
# sig_roohistpdf_met = ROOT.RooHistPdf("sig_roohistpdf_met", "sig_roohistpdf_met", ROOT.RooArgSet(met), sig_roo_template_hist, intOrder=0)


# Signal 1d mll model
mean_mll = ROOT.RooRealVar("mean_mll", "mean_mll", 90, 85, 95)
sigmal_mll = ROOT.RooRealVar("sigmal_mll", "sigmal_mll", 1.7, 0.1, 50)
sigmar_mll = ROOT.RooRealVar("sigmar_mll", "sigmar_mll", 0.5, 0.1, 50)
alphal_mll = ROOT.RooRealVar("alphal_mll","alphal_mll", 0.6, 0.1, 50)
nl_mll = ROOT.RooRealVar("nl_mll", "nl_mll", 9.8, 1, 50)
alphar_mll = ROOT.RooRealVar("alphar_mll","alphar_mll", 0.16, 0.01, 50)
nr_mll = ROOT.RooRealVar("nr_mll", "nr_mll", 73, 1, 200)
sig_dcb_mll = ROOT.RooCrystalBall("sig_dcb_mll", "sig_dcb_mll", mll, mean_mll, sigmal_mll, sigmar_mll, alphal_mll, nl_mll, alphar_mll, nr_mll)

# Signal 2D model: sigtot_mll_met_2dpdf = sig_dcb_mll * pdf_of_spline
sigtot_mll_met_2dpdf = ROOT.RooProdPdf("sigtot_dcb_mll_spline_met", "sigtot_dcb_mll_spline_met", [sig_dcb_mll, pdf_of_spline])

# 2D signal unbinned fit 
sig_result = sigtot_mll_met_2dpdf.fitTo(sigdataset, RF.Save(), SumW2Error=True) #where dataset is RooDataSet
params = sigtot_mll_met_2dpdf.getParameters(sigdataset)

print(params.Print("v"))


###########################################################################
# Declare the background 1D PDFs and multiply them into 2D PDFs
###########################################################################
# Background fake mll component in MET dimension
mu_fakemll_met = ROOT.RooRealVar('mu_fakemll_met', 'mu_fakemll_met', 225, 100, 300)
b_fakemll_met = ROOT.RooRealVar('b_fakemll_met', 'b_fakemll_met', 40, 5, 50)
bkgfakemll_met = ROOT.RooGenericPdf("bkgfakemll_met", "bkgfakemll_met", "1/b_fakemll_met * exp(-(@0 - mu_fakemll_met)/b_fakemll_met - exp(-(@0 - mu_fakemll_met)/b_fakemll_met))",
                        ROOT.RooArgList(met, mu_fakemll_met, b_fakemll_met))
# Background fake mll component in mll dimensio
a_fakemll_mll = ROOT.RooRealVar("a_fakemll_mll", "a_fakemll_mll", -0.03, -1, 1)
bkgfakemll_mll = ROOT.RooExponential("bkgfakemll_mll", "bkgfakemll_mll", mll, a_fakemll_mll)
#Background 2d fakemll model: bkgfakemll_mll_met_2dpdf = bkgfakemll_met * bkgfakemll_mll
bkgfakemll_mll_met_2dpdf = ROOT.RooProdPdf("bkgfakemll_mll_met_2dpdf", "bkgfakemll_mll_met_2dpdf", [bkgfakemll_mll, bkgfakemll_met])


# Background real mll component in MET dimension
mu_realmll_met = ROOT.RooRealVar('mu_realmll_met', 'mu_realmll_met', 50, 30, 400)   # adequate value somewhere around 246 based on hand-drawn plots
b_realmll_met = ROOT.RooRealVar('b_realmll_met', 'b_realmll_met', 40, 20, 100)    # adequate value somewhere around 41 based on hand-drawn plots
bkgrealmll_met = ROOT.RooGenericPdf("bkgrealmll_met", "bkgrealmll_met", "1/b_realmll_met * exp(-(@0 - mu_realmll_met)/b_realmll_met - exp(-(@0 - mu_realmll_met)/b_realmll_met))",
                        ROOT.RooArgList(met, mu_realmll_met, b_realmll_met))
# Background real mll component in mll dimension (parameters taken from initial fit)
zpeak_file = ROOT.TFile.Open("../zpeak_fit/initial_zPeak_fit_result.root", "READ")
zpeak_result = zpeak_file.Get("zPeak_CRZ_fit_result")
zpeak_mean_mll = zpeak_result.floatParsFinal().find("peak_mean_mll")
zpeak_sigmal_mll = zpeak_result.floatParsFinal().find("peak_sigmal_mll")
zpeak_sigmar_mll = zpeak_result.floatParsFinal().find("peak_sigmar_mll")
zpeak_alphal_mll = zpeak_result.floatParsFinal().find("peak_alphal_mll")
zpeak_nl_mll = zpeak_result.floatParsFinal().find("peak_nl_mll")
zpeak_alphar_mll = zpeak_result.floatParsFinal().find("peak_alphar_mll")
zpeak_nr_mll = zpeak_result.floatParsFinal().find("peak_nr_mll")
for v in [zpeak_mean_mll, zpeak_sigmal_mll, zpeak_sigmar_mll, zpeak_alphal_mll, zpeak_nl_mll, zpeak_alphar_mll, zpeak_nr_mll]:
    v.setConstant()
bkgrealmll_mll = ROOT.RooCrystalBall("bkgrealmll_mll", "bkgrealmll_mll", mll, zpeak_mean_mll, zpeak_sigmal_mll, zpeak_sigmar_mll, zpeak_alphal_mll, zpeak_nl_mll, zpeak_alphar_mll, zpeak_nr_mll)

# Background real mll component: 2D PDF (product)
bkgrealmll_mll_met_2dpdf = ROOT.RooProdPdf("bkgrealmll_mll_met_2dpdf", "bkgrealmll_mll_met_2dpdf", [bkgrealmll_mll, bkgrealmll_met])

# Overall 2D bkg model: bkgtot_mll_met_2dpdf = ratio_realmll * bkgrealmll_mll_met_2dpdf + bkgfakemll_mll_met_2dpdf
ratio_realmll = ROOT.RooRealVar("ratio_realmll", "ratio_realmll", 0.1, 0, 1)
bkgtot_mll_met_2dpdf = ROOT.RooAddPdf("bkgtot_mll_met_2dpdf", "bkgtot_mll_met_2dpdf", [bkgrealmll_mll_met_2dpdf, bkgfakemll_mll_met_2dpdf], [ratio_realmll])

# Then fit met projection with ratio_realmll fixed
bkg_result = bkgtot_mll_met_2dpdf.fitTo(bkgdataset, RF.Save(), RF.SumW2Error(True))
params_bkg = bkgtot_mll_met_2dpdf.getParameters(bkgdataset)
print(params_bkg.Print("v"))

w = ROOT.RooWorkspace("workspace", "workspace")
w.Import(sig_roohistpdf_met)
w.Import(pdf_of_spline)

f = ROOT.TFile("fitresult.root", "RECREATE")
sig_result.Write("sig_result")
bkg_result.Write("bkg_result")
zpeak_result.Write("zpeak_result")
w.Write()
f.Close()
