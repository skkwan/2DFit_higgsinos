# test.py
# Diagnostic script: checks why the "total" model curve sits below the
# "background" component curve in the 1D projection plots produced by
# fit2D_bkgOnly_SplusBmodel.py, and fixes the normalization.
#
# Run in ROOT 6.38 (do not do cmsenv)

import ROOT
from ROOT import RooFit as RF
import cmsstyle as CMS
import os

##################################################
# Load pre-fit results (same as original script)
##################################################
resultsfile = ROOT.TFile.Open("../pdf_fit/fitresult.root", "READ")
workspace    = resultsfile.Get("workspace")
sig_results  = resultsfile.Get("sig_result")
bkg_results  = resultsfile.Get("bkg_result")
zpeak_results = resultsfile.Get("zpeak_result")

met = workspace.var("met")
mll = ROOT.RooRealVar("m_ll", "m_ll", 60, 120)

##################################################
# Signal 2D PDF
##################################################
mean_mll   = sig_results.floatParsFinal().find("mean_mll")
sigmal_mll = sig_results.floatParsFinal().find("sigmal_mll")
sigmar_mll = sig_results.floatParsFinal().find("sigmar_mll")
alphal_mll = sig_results.floatParsFinal().find("alphal_mll")
alphar_mll = sig_results.floatParsFinal().find("alphar_mll")
nl_mll     = sig_results.floatParsFinal().find("nl_mll")
nr_mll     = sig_results.floatParsFinal().find("nr_mll")
for v in [mean_mll, sigmal_mll, sigmar_mll, alphal_mll, alphar_mll, nl_mll, nr_mll]:
    v.setConstant(True)

sig_dcb_mll = ROOT.RooCrystalBall("sig_dcb_mll", "sig_dcb_mll",
                                   mll, mean_mll, sigmal_mll, sigmar_mll,
                                   alphal_mll, nl_mll, alphar_mll, nr_mll)

pdf_of_spline = workspace.pdf("pdf_of_spline")

sigtot_mll_met_2dpdf = ROOT.RooProdPdf("sigtot_mll_met_2dpdf", "sigtot_mll_met_2dpdf",
                                        ROOT.RooArgList(sig_dcb_mll, pdf_of_spline))

##################################################
# Background 2D PDF
##################################################
mu_fakemll_met = bkg_results.floatParsFinal().find("mu_fakemll_met")
b_fakemll_met  = bkg_results.floatParsFinal().find("b_fakemll_met")
a_fakemll_mll  = bkg_results.floatParsFinal().find("a_fakemll_mll")
for v in [mu_fakemll_met, b_fakemll_met, a_fakemll_mll]:
    v.setConstant(True)

bkgfakemll_met = ROOT.RooGenericPdf(
    "bkgfakemll_met", "bkgfakemll_met",
    "1/b_fakemll_met * exp(-(@0 - mu_fakemll_met)/b_fakemll_met"
    " - exp(-(@0 - mu_fakemll_met)/b_fakemll_met))",
    ROOT.RooArgList(met, mu_fakemll_met, b_fakemll_met))

bkgfakemll_mll = ROOT.RooExponential("bkgfakemll_mll", "bkgfakemll_mll",
                                      mll, a_fakemll_mll)

bkgfakemll_mll_met_2dpdf = ROOT.RooProdPdf("bkgfakemll_mll_met_2dpdf",
                                             "bkgfakemll_mll_met_2dpdf",
                                             ROOT.RooArgList(bkgfakemll_mll, bkgfakemll_met))

mu_realmll_met = bkg_results.floatParsFinal().find("mu_realmll_met")
b_realmll_met  = bkg_results.floatParsFinal().find("b_realmll_met")
for v in [mu_realmll_met, b_realmll_met]:
    v.setConstant(True)

bkgrealmll_met = ROOT.RooGenericPdf(
    "bkgrealmll_met", "bkgrealmll_met",
    "1/b_realmll_met * exp(-(@0 - mu_realmll_met)/b_realmll_met"
    " - exp(-(@0 - mu_realmll_met)/b_realmll_met))",
    ROOT.RooArgList(met, mu_realmll_met, b_realmll_met))

zpeak_mean_mll   = zpeak_results.floatParsFinal().find("peak_mean_mll")
zpeak_sigmal_mll = zpeak_results.floatParsFinal().find("peak_sigmal_mll")
zpeak_sigmar_mll = zpeak_results.floatParsFinal().find("peak_sigmar_mll")
zpeak_alphal_mll = zpeak_results.floatParsFinal().find("peak_alphal_mll")
zpeak_nl_mll     = zpeak_results.floatParsFinal().find("peak_nl_mll")
zpeak_alphar_mll = zpeak_results.floatParsFinal().find("peak_alphar_mll")
zpeak_nr_mll     = zpeak_results.floatParsFinal().find("peak_nr_mll")
for v in [zpeak_mean_mll, zpeak_sigmal_mll, zpeak_sigmar_mll,
          zpeak_alphal_mll, zpeak_nl_mll, zpeak_alphar_mll, zpeak_nr_mll]:
    v.setConstant(True)

bkgrealmll_mll = ROOT.RooCrystalBall("bkgrealmll_mll", "bkgrealmll_mll",
                                      mll, zpeak_mean_mll, zpeak_sigmal_mll, zpeak_sigmar_mll,
                                      zpeak_alphal_mll, zpeak_nl_mll, zpeak_alphar_mll, zpeak_nr_mll)

bkgrealmll_mll_met_2dpdf = ROOT.RooProdPdf("bkgrealmll_mll_met_2dpdf",
                                             "bkgrealmll_mll_met_2dpdf",
                                             ROOT.RooArgList(bkgrealmll_mll, bkgrealmll_met))

ratio_realmll = ROOT.RooRealVar("ratio_realmll", "ratio_realmll", 0.1, 0, 1)

bkgtot_mll_met_2dpdf = ROOT.RooAddPdf("bkgtot_mll_met_2dpdf", "bkgtot_mll_met_2dpdf",
                                       ROOT.RooArgList(bkgrealmll_mll_met_2dpdf,
                                                       bkgfakemll_mll_met_2dpdf),
                                       ratio_realmll)

##################################################
# Dataset
##################################################
combfile = ROOT.TFile.Open("backgrounds_for_2D_fit.root", "READ")
combtree = combfile.Get("event_tree")

weightVar = ROOT.RooRealVar("weight_nominal", "weight_nominal", -1, 1)
variables = ROOT.RooArgSet(mll, met, weightVar)
combdataset = ROOT.RooDataSet("combdataset", "Signal + background MC",
                               variables,
                               RF.Import(combtree),
                               RF.WeightVar(weightVar))

n_total = combdataset.sumEntries()
print(f"\nDataset: {combdataset.numEntries()} entries, sum of weights = {n_total:.2f}")

##################################################
# Extended model + fit
##################################################
n_sig = ROOT.RooRealVar("n_sig", "Signal yield",     n_total * 0.1, 0, n_total * 10)
n_bkg = ROOT.RooRealVar("n_bkg", "Background yield", n_total * 0.9, 0, n_total * 10)

model_2dpdf = ROOT.RooAddPdf("model_2dpdf", "Signal + Background 2D model",
                              ROOT.RooArgList(sigtot_mll_met_2dpdf, bkgtot_mll_met_2dpdf),
                              ROOT.RooArgList(n_sig, n_bkg))

fit_result = model_2dpdf.fitTo(combdataset,
                                RF.Extended(True),
                                RF.Save(True),
                                RF.SumW2Error(True),
                                RF.PrintLevel(1))

##################################################
# CHECK 1: Yield diagnostics
# If n_sig + n_bkg < n_bkg the total curve will be below the background
# curve. This happens when n_sig is negative (or the ratio_realmll
# normalization inflates the component relative to the full model).
##################################################
print("\n=== Yield check ===")
print(f"  n_sig            = {n_sig.getVal():.4f}  (lower bound = {n_sig.getMin():.1f})")
print(f"  n_bkg            = {n_bkg.getVal():.4f}")
print(f"  n_sig + n_bkg    = {n_sig.getVal() + n_bkg.getVal():.4f}  <- total curve normalisation")
print(f"  n_total (data)   = {n_total:.4f}")

if n_sig.getVal() < 0:
    print("  [ISSUE] n_sig is NEGATIVE: total = sig+bkg < bkg, so blue sits below green.")
elif abs(n_sig.getVal() - n_sig.getMin()) < 1e-3:
    print("  [ISSUE] n_sig is at its lower bound — fitter tried to go negative but was clipped.")
else:
    print("  n_sig is positive and not at boundary.")

##################################################
# CHECK 2: Component normalization via explicit integral comparison
# RooFit's RF.Components plots the sub-PDF scaled to the *full* model norm,
# but nested RooAddPdfs can be renormalized relative to their own internal
# coefficient, inflating or deflating the component curve.
##################################################
print("\n=== Normalization check (1D integrals over mll) ===")
mll_obs = ROOT.RooArgSet(mll)

# Integral of the full model PDF over mll (should be 1 for a normalised PDF)
full_integral = model_2dpdf.createIntegral(mll_obs)
bkg_integral  = bkgtot_mll_met_2dpdf.createIntegral(mll_obs)
print(f"  Integral of model_2dpdf over mll        = {full_integral.getVal():.6f}")
print(f"  Integral of bkgtot_mll_met_2dpdf over mll = {bkg_integral.getVal():.6f}")
print("  (Both should be 1.0 for properly normalised unit PDFs; the visual")
print("   height is then set by the event yield scale factor.)")

print("\n=== Fit summary ===")
print(f"  Fit status: {fit_result.status()}  (0 = OK)")
print(f"  EDM:        {fit_result.edm():.3e}")
fit_result.Print("v")

##################################################
# FIX: Projection plots with explicit RF.Normalization on each component
# so the vertical scale is exactly n_sig / n_bkg events, regardless of
# how RooFit internally re-normalises nested RooAddPdf coefficients.
##################################################
ROOT.gROOT.SetBatch(True)
CMS.setCMSStyle()
CMS.SetExtraText("Simulation Preliminary")
CMS.cms_lumi = "(13 TeV)"
CMS.cms_energy = ""
CMS.lumi_13TeV = ""

obs_configs = [
    (mll, 30,   60,  120, "m(ll) [GeV]", "test_proj_mll"),
    (met, 40,    0, 1200, "MET [GeV]",   "test_proj_met"),
]

n_sig_val  = n_sig.getVal()
n_bkg_val  = n_bkg.getVal()
n_tot_val  = n_sig_val + n_bkg_val

for obs, nBins, xmin, xmax, xlabel, plotname in obs_configs:
    frame = obs.frame(RF.Bins(nBins), RF.Range(xmin, xmax), RF.Title(""))

    # Total: normalised to n_sig + n_bkg
    model_2dpdf.plotOn(frame,
                       RF.Name("total"),
                       RF.Normalization(n_tot_val, ROOT.RooAbsReal.NumEvent),
                       RF.LineColor(ROOT.kBlue),
                       RF.LineWidth(2))

    # Signal component: normalised to exactly n_sig events
    model_2dpdf.plotOn(frame,
                       RF.Components("sigtot_mll_met_2dpdf"),
                       RF.Name("signal"),
                       RF.Normalization(n_sig_val, ROOT.RooAbsReal.NumEvent),
                       RF.LineColor(ROOT.kRed),
                       RF.LineStyle(ROOT.kDashed),
                       RF.LineWidth(2))

    # Background component: normalised to exactly n_bkg events
    model_2dpdf.plotOn(frame,
                       RF.Components("bkgtot_mll_met_2dpdf"),
                       RF.Name("background"),
                       RF.Normalization(n_bkg_val, ROOT.RooAbsReal.NumEvent),
                       RF.LineColor(ROOT.kGreen + 2),
                       RF.LineStyle(ROOT.kDashed),
                       RF.LineWidth(2))

    # Data
    combdataset.plotOn(frame,
                       RF.Binning(nBins),
                       RF.Name("data"),
                       RF.MarkerColor(ROOT.kBlack),
                       RF.LineColor(ROOT.kBlack),
                       RF.MarkerStyle(ROOT.kFullCircle),
                       RF.MarkerSize(0.8))

    pull_hist = frame.pullHist("data", "total")
    pull_hist.SetMarkerStyle(ROOT.kFullCircle)
    pull_hist.SetMarkerSize(0.8)

    frame_pull = obs.frame(RF.Bins(nBins), RF.Range(xmin, xmax), RF.Title(""))
    frame_pull.addPlotable(pull_hist, "P")

    c = ROOT.TCanvas(plotname, plotname, 800, 800)
    pad_main = ROOT.TPad("pad_main", "", 0, 0.28, 1, 1.0)
    pad_pull = ROOT.TPad("pad_pull", "", 0, 0.00, 1, 0.28)
    pad_main.SetBottomMargin(0.02)
    pad_main.SetTopMargin(0.07)
    pad_main.SetLeftMargin(0.12)
    pad_pull.SetTopMargin(0.04)
    pad_pull.SetBottomMargin(0.35)
    pad_pull.SetLeftMargin(0.12)
    pad_main.Draw()
    pad_pull.Draw()

    pad_main.cd()
    frame.GetXaxis().SetLabelSize(0)
    frame.GetXaxis().SetTitleSize(0)
    frame.GetYaxis().SetTitle("Events / bin")
    frame.GetYaxis().SetTitleSize(0.055)
    frame.GetYaxis().SetTitleOffset(1.0)
    frame.Draw()

    leg = ROOT.TLegend(0.55, 0.55, 0.90, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.042)
    leg.AddEntry(frame.findObject("data"),      "MC (signal + bkg)", "PE")
    leg.AddEntry(frame.findObject("total"),      "Total model", "L")
    leg.AddEntry(frame.findObject("signal"),     f"Signal  (n_{{sig}} = {n_sig_val:.1f})", "L")
    leg.AddEntry(frame.findObject("background"), f"Background  (n_{{bkg}} = {n_bkg_val:.1f})", "L")
    leg.Draw()

    CMS.CMS_lumi(pad_main, iPosX=0)

    pad_pull.cd()
    frame_pull.GetYaxis().SetTitle("Pull")
    frame_pull.GetYaxis().SetNdivisions(5)
    frame_pull.GetYaxis().SetTitleSize(0.13)
    frame_pull.GetYaxis().SetTitleOffset(0.38)
    frame_pull.GetYaxis().SetLabelSize(0.10)
    frame_pull.GetXaxis().SetTitle(xlabel)
    frame_pull.GetXaxis().SetTitleSize(0.14)
    frame_pull.GetXaxis().SetTitleOffset(0.85)
    frame_pull.GetXaxis().SetLabelSize(0.11)
    frame_pull.SetMaximum(4.9)
    frame_pull.SetMinimum(-4.9)
    frame_pull.Draw()

    zero_line = ROOT.TLine(xmin, 0, xmax, 0)
    zero_line.SetLineColor(ROOT.kBlue)
    zero_line.SetLineWidth(1)
    zero_line.Draw("SAME")

    c.SaveAs(f"{plotname}.pdf")
    c.SaveAs(f"{plotname}.png")
    print(f"Created {plotname}.pdf / .png")
