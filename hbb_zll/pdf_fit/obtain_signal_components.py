# Run in ROOT 6.38 (do not do cmsenv)

import os
import ROOT
from ROOT import RooFit as RF
import numpy as np
from array import array
import cmsstyle as CMS

doLog = True


def plotSignalFit(name, rooVar, dataset, pdf, dataLabel, fitLabel, plotname, nFloatParams=2, outdir="", getOverflow=True, doLog=True):
    nBins = 50
    xmin = rooVar.getMin()
    xmax = rooVar.getMax()

    scale_suffix = "_log" if doLog else "_lin"
    frame = rooVar.frame(nBins)

    leg = CMS.cmsLeg(0.3, 0.89 - 0.05 * 4, 0.9, 0.89, textSize=0.04)
    CMS.SetLumi("")

    data_hist = dataset.createHistogram("histo_" + plotname + scale_suffix, rooVar, ROOT.RooFit.Binning(nBins, xmin, xmax))
    if getOverflow:
        data_hist.SetBinContent(data_hist.GetNbinsX(), data_hist.GetBinContent(data_hist.GetNbinsX()) + data_hist.GetBinContent(data_hist.GetNbinsX() + 1))
    pdf_hist = pdf.createHistogram("hpdf_" + plotname + scale_suffix, rooVar, ROOT.RooFit.Binning(nBins))
    if getOverflow:
        pdf_hist.SetBinContent(pdf_hist.GetNbinsX(), pdf_hist.GetBinContent(pdf_hist.GetNbinsX()) + pdf_hist.GetBinContent(pdf_hist.GetNbinsX() + 1))
    y_min = 0
    y_max = 1.8 * max(data_hist.GetMaximum(), pdf_hist.GetMaximum())
    if doLog:
        y_min = 1e-10
        y_max = y_max * 1000

    canv = CMS.cmsDiCanvas("canv_" + plotname + scale_suffix, x_min=xmin, x_max=xmax, y_min=y_min, y_max=y_max,
                           r_min=0, r_max=2,
                           nameXaxis=f"{name} / GeV",
                           nameYaxis="Shape (A.U.)",
                           nameRatio="MC/Pred",
                           square=CMS.kSquare, iPos=0)
    canv.SetRightMargin(0.05)
    CMS.UpdatePad(canv)

    canv.cd(1)
    if doLog:
        ROOT.gPad.SetLogy()
    CMS.UpdatePad(canv)

    dataset.plotOn(frame, ROOT.RooFit.Name("data_" + plotname),
                   ROOT.RooFit.LineColor(ROOT.TColor.GetColor("#5790fc")),
                   ROOT.RooFit.LineWidth(2),
                   ROOT.RooFit.MarkerColor(ROOT.TColor.GetColor("#5790fc")),
                   ROOT.RooFit.MarkerSize(1),
                   ROOT.RooFit.Binning(nBins))
    leg.AddEntry(frame.findObject("data_" + plotname), dataLabel)

    pdf.plotOn(frame, ROOT.RooFit.Name("pdf_" + plotname),
               ROOT.RooFit.LineColor(ROOT.TColor.GetColor("#f89c20")),
               ROOT.RooFit.LineWidth(2),
               ROOT.RooFit.LineStyle(1),
               ROOT.RooFit.MarkerSize(0),
               ROOT.RooFit.Binning(nBins))
    leg.AddEntry(frame.findObject("pdf_" + plotname), fitLabel)

    roo_curve = frame.getCurve("pdf_" + plotname)
    chi2_per_ndf = frame.chiSquare("pdf_" + plotname, "data_" + plotname, nFloatParams)
    leg.SetHeader(f"2018 SR: (650, 1) GeV signal (#chi^{{2}}/ndf = {chi2_per_ndf:.2f})")
    frame.Draw("SAME")

    # Ratio plot
    canv.cd(2)
    data_ratio = data_hist.Clone()
    prediction = data_hist.Clone()
    for i in range(1, data_ratio.GetNbinsX() + 1):
        thisXval = data_ratio.GetBinCenter(i)
        pdfY = roo_curve.Eval(thisXval)
        prediction.SetBinContent(i, pdfY)
    data_ratio.Divide(prediction)
    data_ratio.SetMarkerColor(ROOT.kBlack)
    data_ratio.SetLineColor(ROOT.kBlack)
    CMS.cmsObjectDraw(data_ratio, "E", MarkerStyle=ROOT.kFullCircle)
    unitLine = ROOT.TLine(xmin, 1.0, xmax, 1.0)
    unitLine.SetLineColor(ROOT.kBlack)
    unitLine.SetLineWidth(1)
    unitLine.Draw("SAME")

    canv.cd(1)
    CMS.cmsObjectDraw(leg)
    CMS.UpdatePad(canv)

    fname = plotname + ("-log" if doLog else "")
    canv.SaveAs(f"{fname}.pdf")
    canv.SaveAs(f"{fname}.png")
    if outdir:
        os.system(f"mv {fname}.* {outdir}")

    del canv
    del leg

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
met = ROOT.RooRealVar("met", "met", 200, 1200)

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
sig_met_hist = ROOT.TH1D("sig_met_hist", "sig_met_hist", 50, 200, 1200)
sigtree.Draw("met >> sig_met_hist")
# Convert TH1 into RooDataHist
sig_roo_template_hist = ROOT.RooDataHist("sig_roo_template_hist", "sig_roo_template_hist", met, sig_met_hist)
# Create a RooHistPdf based on the RooFit histogram
sig_roohistpdf_met = ROOT.RooHistPdf("sig_roohistpdf_met", "sig_roohistpdf_met", met, sig_roo_template_hist, intOrder=0)

# Build spline-based signal met PDF from signal met histogram
nBins = sig_met_hist.GetNbinsX()
centers = np.array([sig_met_hist.GetBinCenter(i) for i in range(1, nBins+1)], dtype=np.float64)
vals    = np.array([max(float(sig_met_hist.GetBinContent(i)), 0.0) for i in range(1, nBins+1)], dtype=np.float64)

knot_avg_halfwidth_bins = 2
min_y = 1e-8
tail_thresh=1e-3

# Restrict knots to non-zero bins so spline is never built over the MET values where the signal is less than min_y
mask = np.array(vals > min_y)
centers_nz = centers[mask]
vals_nz    = vals[mask]
first_knot_x = float(centers_nz[0])
print(f"First non-zero bin: {first_knot_x:.1f} GeV  (masked out {np.sum(~mask)} zero bins)")

knot_x = make_knot_x(
    x_min=centers[0],
    x_max=centers[-1],
    n_knots=240,
    power=2, #higher number = more dense at low mass
    centers=centers,
)
print("knot_x:", knot_x)
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


# knot_x, knot_y = prune_knots_in_tiny_tail(knot_x, knot_y, tail_thresh=tail_thresh)
# knot_y_use = np.maximum(knot_y, min_y).astype(np.double)
vx = ROOT.std.vector('double')(knot_x)
vy_use = ROOT.std.vector('double')(knot_y)

print("Zipped:", list(zip(vx, vy_use)))

spline = ROOT.RooSpline(
        "spline", "spline",
        met,
        vx,
        vy_use,
        order=3,
    )

# pdf_of_spline = ROOT.RooGenericPdf(
#     "pdf_of_spline", "pdf_of_spline",
#     "max(@0, 1e-12)",    # prevent the interpolation from going 0 or negative
#     ROOT.RooArgList(spline)
# )

# Step clamp: PDF is exactly 1e-12 below the first knot, preventing wild extrapolation
pdf_of_spline = ROOT.RooGenericPdf(
    "pdf_of_spline", "pdf_of_spline",
    f"(@1 < {first_knot_x}) ? 1e-12 : max(@0, 1e-12)",
    ROOT.RooArgList(spline, met)
)


# # Also keep a RooHistPdf for reference
# sig_roo_template_hist = ROOT.RooDataHist("sig_roo_template_hist", "sig_roo_template_hist", ROOT.RooArgSet(met), sig_met_hist)
# sig_roohistpdf_met = ROOT.RooHistPdf("sig_roohistpdf_met", "sig_roohistpdf_met", ROOT.RooArgSet(met), sig_roo_template_hist, intOrder=0)

# Signal 1d mll model
mean_mll = ROOT.RooRealVar("mean_mll", "mean_mll", 90, 85, 95)
sigmal_mll = ROOT.RooRealVar("sigmal_mll", "sigmal_mll", 1.7, 0.1, 50)
sigmar_mll = ROOT.RooRealVar("sigmar_mll", "sigmar_mll", 0.5, 0.1, 50)
alphal_mll = ROOT.RooRealVar("alphal_mll","alphal_mll", 2.4, 1, 50)
nl_mll = ROOT.RooRealVar("nl_mll", "nl_mll", 2.8, 0.5, 20)
alphar_mll = ROOT.RooRealVar("alphar_mll","alphar_mll", 2.4, 1, 50)
nr_mll = ROOT.RooRealVar("nr_mll", "nr_mll", 2.25, 0.5, 20)
sig_dcb_mll = ROOT.RooCrystalBall("sig_dcb_mll", "sig_dcb_mll", mll, mean_mll, sigmal_mll, sigmar_mll, alphal_mll, nl_mll, alphar_mll, nr_mll)

# Signal 2D model: sigtot_mll_met_2dpdf = sig_dcb_mll * pdf_of_spline
sigtot_mll_met_2dpdf = ROOT.RooProdPdf("sigtot_dcb_mll_spline_met", "sigtot_dcb_mll_spline_met", [sig_dcb_mll, pdf_of_spline])

# 2D signal unbinned fit 
sig_result = sigtot_mll_met_2dpdf.fitTo(sigdataset, RF.Save(), SumW2Error=True) #where dataset is RooDataSet
params = sigtot_mll_met_2dpdf.getParameters(sigdataset)

print(params.Print("v"))


w = ROOT.RooWorkspace("workspace", "workspace")
w.Import(sig_roohistpdf_met)
w.Import(pdf_of_spline)
# w.Import(spline)

f = ROOT.TFile("fitresult_signal.root", "RECREATE")
sig_result.Write("sig_result")
# bkg_result.Write("bkg_result")
# zpeak_result.Write("zpeak_result")
w.Write()
f.Close()

for doLog in [True, False]:
    plotSignalFit("MET", met, sigdataset, pdf_of_spline,
                  "Signal MC",
                  "Spline fit",
                  "sig_met_spline",
                  doLog=doLog)
    plotSignalFit("m(ll)", mll, sigdataset, sig_dcb_mll,
                  "Signal MC",
                  f"DCB fit (#mu={mean_mll.getVal():.1f}#pm{mean_mll.getError():.1f}, #sigma_L={sigmal_mll.getVal():.1f}#pm{sigmal_mll.getError():.1f})",
                  "sig_mll_dcb",
                  nFloatParams=7,
                  doLog=doLog)

os.system("mv sig*.pdf /eos/user/s/skkwan/www/higgsino/studies/mll-MET-fit-2D/signal_shapes")
os.system("mv sig*.png /eos/user/s/skkwan/www/higgsino/studies/mll-MET-fit-2D/signal_shapes")
