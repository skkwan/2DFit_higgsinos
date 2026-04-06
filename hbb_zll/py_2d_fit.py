# Run in ROOT 6.38 (do not do cmsenv)

import ROOT
from ROOT import RooFit as RF
import numpy as np
from array import array

# From https://gitlab.cern.ch/cms-l1-ad/coffea-dask-axol1tl-studies/-/blob/master/prepareDatacards.py?ref_type=heads 
def prune_knots_in_tiny_tail(knot_x, knot_y, tail_thresh=1e-4, min_keep=10):
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


##### DEFINE FIT OBSERVABLES ####
mll = ROOT.RooRealVar("m_ll", "m_ll", 60, 120) 
met = ROOT.RooRealVar("met", "met", 0, 1200)


###### Retrieve signal dataset from signal root file 
sigfilepath = 'snapshot_TChiZH_650_1_cat_0_batch_0_channel_mm_SR_mll_MET_fit_scheme.root'
sigfile = ROOT.TFile.Open(sigfilepath, "READ")
sigtree = sigfile.Get("event_tree")
weightXyear = ROOT.RooRealVar("weight_nominal_mm", "weight_nominal_mm", -1, 1)
variables = ROOT.RooArgSet(mll, met, weightXyear)
sigdataset = ROOT.RooDataSet("sigdataset", "sigdataset", variables, ROOT.RooFit.Import(sigtree), ROOT.RooFit.WeightVar(weightXyear))


# # ##### SIGNAL FIT to signal MC file ######

#Signal 1d met model. Formerly sigmoid in Meraj's notebook, we try a DCB and also a Gamma
a_met = ROOT.RooRealVar('a_met', 'a_met', 10, 0, 3000)
b_met = ROOT.RooRealVar('b_met', 'b_met', 10, 0, 1000)  
c_met = ROOT.RooRealVar('c_met', 'c_met', 0.5, 0, 1000)  
e_met = ROOT.RooRealVar('e_met', 'e_met', 1, 0.2, 100) 
sig_smoid_met = ROOT.RooGenericPdf('sig_smoid_met', '(1-exp(-c_met*met))/(1 + exp((met^e_met-a_met)/b_met))', ROOT.RooArgList(met, a_met, b_met, c_met, e_met))

# Get TH1 of signal met: https://github.com/guitargeek/hasco-2023-root/blob/main/notebooks/roofit-tutorial-01.ipynb 
sig_met_hist = ROOT.TH1D("sig_met_hist", "sig_met_hist", 120, 0, 1200)
sigtree.Draw("met >> sig_met_hist")
# Convert TH1 into RooDataHist
sig_roo_template_hist = ROOT.RooDataHist("sig_roo_template_hist", "sig_roo_template_hist", met, sig_met_hist)
# Create a RooHistPdf based on the RooFit histogram
sig_roohistpdf_met = ROOT.RooHistPdf("sig_roohistpdf_met", "sig_roohistpdf_met", met, sig_roo_template_hist, intOrder=0)


hIn = sig_met_hist 
nBins = hIn.GetNbinsX()
centers = np.array([hIn.GetBinCenter(i) for i in range(1, nBins+1)], dtype=np.float64)
vals    = np.array([max(float(hIn.GetBinContent(i)), 0.0) for i in range(1, nBins+1)], dtype=np.float64)

knot_avg_halfwidth_bins = 4
min_y = 1e-4
l=1.0
tail_thresh=1e-4
knot_x = make_knot_x(
    x_min=centers[0],
    x_max=centers[-1],
    n_knots=60,
    power = 10.0, #higher number = more dense at low mass
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
vx = ROOT.std.vector('double')(knot_x)
vy = ROOT.std.vector('double')(knot_y)

# knot_x, knot_y = prune_knots_in_tiny_tail(knot_x, knot_y, tail_thresh=tail_thresh, min_keep=2)
knot_y_use = np.maximum(knot_y, min_y).astype(np.double)
vy_use = ROOT.std.vector('double')(knot_y_use)

spline = ROOT.RooSpline(
        "spline", "spline",
        met,
        vx,
        vy_use,
        order=3
    )

pdf_of_spline = ROOT.RooGenericPdf(
    "pdf_of_spline", "pdf_of_spline",
    "max(@0, 1e-9)",    # prevent the interpolation from going 0 or negative 
    ROOT.RooArgList(spline)
)

# Load the TH1 into here


# TODO: testing DCB for signal met
mean_met = ROOT.RooRealVar("mean_met", "mean_met", 400, 300, 500)
sigmal_met = ROOT.RooRealVar("sigmal_met", "sigmal_met", 80, 30, 200)
sigmar_met = ROOT.RooRealVar("sigmar_met", "sigmar_met", 2, 0.5, 10)
alphal_met = ROOT.RooRealVar("alphal_met","alphal_met", 4, 0.01, 10)
nl_met = ROOT.RooRealVar("nl_met", "nl_met", 3, 1, 10)
alphar_met = ROOT.RooRealVar("alphar_met","alphar_met", 5, 0.01, 10)
nr_met = ROOT.RooRealVar("nr_met", "nr_met", 3, 1, 10)
sig_dcb_met = ROOT.RooCrystalBall("sig_dcb_met", "sig_dcb_met", met, mean_met, sigmal_met, sigmar_met, alphal_met, nl_met, alphar_met, nr_met)

# https://root.cern.ch/doc/v638/classRooGamma.html
gamma_met = ROOT.RooRealVar("gamma_met", "gamma_met", 200, 10, 500) # gamma in ROOT = alpha on wikipedia
beta_met = ROOT.RooRealVar("beta_met", "beta_met", 1.0, 0.5, 2.5) # wikipedia: beta = 0.5 to 1.0, beta in wikipedia = (1/beta) in ROOT, so 
mu_met = ROOT.RooRealVar("mu_met", "mu_met", 0, 0, 1)
sig_gamma_met = ROOT.RooGamma("sig_gamma_met", "sig_gamma_met", met, gamma_met, beta_met, mu_met)

#Signal 1d mll model
mean_mll = ROOT.RooRealVar("mean_mll", "mean_mll", 90, 85, 95)
sigmal_mll = ROOT.RooRealVar("sigmal_mll", "sigmal_mll", 1.7, 0.1, 50)
sigmar_mll = ROOT.RooRealVar("sigmar_mll", "sigmar_mll", 0.5, 0.1, 50)
alphal_mll = ROOT.RooRealVar("alphal_mll","alphal_mll", 0.6, 0.1, 50)
nl_mll = ROOT.RooRealVar("nl_mll", "nl_mll", 9.8, 1, 50)
alphar_mll = ROOT.RooRealVar("alphar_mll","alphar_mll", 0.16, 0.01, 50)
nr_mll = ROOT.RooRealVar("nr_mll", "nr_mll", 73, 1, 200)
sig_dcb_mll = ROOT.RooCrystalBall("sig_dcb_mll", "sig_dcb_mll", mll, mean_mll, sigmal_mll, sigmar_mll, alphal_mll, nl_mll, alphar_mll, nr_mll)

#Signal 2D model: sigtot_mll_met_2dpdf = sig_smoid_met * spline
sigtot_mll_met_2dpdf = ROOT.RooProdPdf("sigtot_dcb_mll_spline_met", "sigtot_dcb_mll_spline_met", [sig_dcb_mll, pdf_of_spline])

###### 2D signal fit 
signal_result = sigtot_mll_met_2dpdf.fitTo(sigdataset, RF.Save(), SumW2Error=True) #where dataset is RooDataSet
params = signal_result.floatParsFinal()

print(params)

#### BACKGROUND FIT to data in control region ######

# Background fake mll model in MET dimension
mu_fakemll_met = ROOT.RooRealVar('mu_fakemll_met', 'mu_fakemll_met', 225, 100, 300) 
b_fakemll_met = ROOT.RooRealVar('b_fakemll_met', 'b_fakemll_met', 40, 5, 50) 
bkgfakemll_met = ROOT.RooGenericPdf("bkgfakemll_met", "bkgfakemll_met", "1/b_fakemll_met * exp(-(@0 - mu_fakemll_met)/b_fakemll_met - exp(-(@0 - mu_fakemll_met)/b_fakemll_met))",
                        ROOT.RooArgList(met, mu_fakemll_met, b_fakemll_met))  
# Background fake mll model in mll dimension: falling exponential (using this for now)
a_fakemll_mll = ROOT.RooRealVar("a_fakemll_mll", "a_fakemll_mll", -0.03, -1, 1) 
bkgfakemll_mll = ROOT.RooExponential("bkgfakemll_mll", "bkgfakemll_mll", mll, a_fakemll_mll)
#Background 2d fakemll model: bkgfakemll_mll_met_2dpdf = bkgfakemll_met * bkgfakemll_mll
bkgfakemll_mll_met_2dpdf = ROOT.RooProdPdf("bkgfakemll_mll_met_2dpdf", "bkgfakemll_mll_met_2dpdf", [bkgfakemll_mll, bkgfakemll_met])


#Background real mll model in met dimension
mu_realmll_met = ROOT.RooRealVar('mu_realmll_met', 'mu_realmll_met', 50, 30, 400)   # adequate value somewhere around 246 based on hand-drawn plots
b_realmll_met = ROOT.RooRealVar('b_realmll_met', 'b_realmll_met', 40, 20, 100)    # adequate value somewhere around 41 based on hand-drawn plots
bkgrealmll_met = ROOT.RooGenericPdf("bkgrealmll_met", "bkgrealmll_met", "1/b_realmll_met * exp(-(@0 - mu_realmll_met)/b_realmll_met - exp(-(@0 - mu_realmll_met)/b_realmll_met))",
                        ROOT.RooArgList(met, mu_realmll_met, b_realmll_met))  
# Background real mll model in mll dimension: use a simple Gaussian
bkg_mean_mll = ROOT.RooRealVar("bkg_mean_mll", "bkg_mean_mll", 90, 85, 95)
bkg_sigma_mll = ROOT.RooRealVar("bkg_sigma_mll", "bkg_sigma_mll", 2, 0.01, 10)
bkgrealmll_sigma_mll = ROOT.RooGaussian("bkg_gaus_mll", "bkg_gaus_mll", mll, bkg_mean_mll, bkg_sigma_mll)

#Background 2d realmll model: bkgrealmll_mll_met_2dpdf = bkgrealmll_met * bkgrealmll_dcb_mll
bkgrealmll_mll_met_2dpdf = ROOT.RooProdPdf("bkgrealmll_mll_met_2dpdf", "bkgrealmll_mll_met_2dpdf", [bkgrealmll_sigma_mll, bkgrealmll_met])


#Overall 2D bkg model: bkgtot_mll_met_2dpdf = bkgfakemet_mll_met_2dpdf + ratio_realmll * bkgrealmll_mll_met_2dpdf
ratio_realmll = ROOT.RooRealVar("ratio_realmll", "ratio_realmll", 0.1, 0, 1)
bkgtot_mll_met_2dpdf = ROOT.RooAddPdf("bkgtot_mll_met_2dpdf", "bkgtot_mll_met_2dpdf", [bkgrealmll_mll_met_2dpdf, bkgfakemll_mll_met_2dpdf], [ratio_realmll])

###### Retrive cr data root file ########
crfilepath = 'backgrounds.root'
crfile = ROOT.TFile.Open(crfilepath, "READ")
crtree = crfile.Get("event_tree")
variables = ROOT.RooArgSet(mll, met)
weight = ROOT.RooRealVar("weight_nominal_mm", "weight_nominal_mm", -1, 1)
crdataset = ROOT.RooDataSet("crdataset", "crdataset", variables, ROOT.RooFit.Import(crtree), ROOT.RooFit.WeightVar(weight))

### B only 2D fit to cr root file
bkg_result = bkgtot_mll_met_2dpdf.fitTo(crdataset, RF.Save(), SumW2Error=True) #where dataset is RooDataSet
params = bkg_result.floatParsFinal()
print(params)

w = ROOT.RooWorkspace("workspace", "workspace")
w.Import(sig_roohistpdf_met)
w.Import(pdf_of_spline)

f = ROOT.TFile("fitresult.root", "RECREATE")
signal_result.Write("signal_result")
bkg_result.Write("bkg_result")
w.Write()
f.Close()

