# fit2D_signal_plus_background.py
# Run in ROOT 6.38 (do not do cmsenv)
#
# Performs a 2D (mll, MET) extended maximum-likelihood fit to
# signal_plus_background_total_MC.root.
#
# All shape parameters are fixed from fitresult_modified.root.
# The only floating parameters are n_sig and n_bkg.

import ROOT
from ROOT import RooFit as RF
import cmsstyle as CMS
import os

##################################################
# Load pre-fit results
##################################################
resultsfile = ROOT.TFile.Open("../pdf_fit/fitresult.root", "READ")
workspace    = resultsfile.Get("workspace")
sig_results  = resultsfile.Get("sig_result")
bkg_results  = resultsfile.Get("bkg_result")
zpeak_results = resultsfile.Get("zpeak_result")

##################################################
# Observables
# Use the workspace's met variable so that pdf_of_spline
# (which was built against it) shares the same object.
##################################################
met = workspace.var("met")  
mll = ROOT.RooRealVar("m_ll", "m_ll", 60, 120)

##################################################
# Signal 2D PDF  (fixed shape)
# sigtot_mll_met_2dpdf = sig_dcb_mll (mll) x pdf_of_spline (MET)
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

# Spline-based MET PDF retrieved from the workspace (shape fixed by construction)
pdf_of_spline = workspace.pdf("pdf_of_spline")

sigtot_mll_met_2dpdf = ROOT.RooProdPdf("sigtot_mll_met_2dpdf", "sigtot_mll_met_2dpdf",
                                        ROOT.RooArgList(sig_dcb_mll, pdf_of_spline))

##################################################
# Background 2D PDF  (fixed shape)
# bkgtot_mll_met_2dpdf = ratio_realmll * bkgrealmll_2dpdf
#                      + (1 - ratio_realmll) * bkgfakemll_2dpdf
##################################################

# --- fake-mll component ---
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

# --- real-mll component (Z peak) ---
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

# --- mixture (ratio_realmll is FLOATING) ---
ratio_realmll = ROOT.RooRealVar("ratio_realmll", "ratio_realmll", 0.1, 0, 1)

bkgtot_mll_met_2dpdf = ROOT.RooAddPdf("bkgtot_mll_met_2dpdf", "bkgtot_mll_met_2dpdf",
                                       ROOT.RooArgList(bkgrealmll_mll_met_2dpdf,
                                                       bkgfakemll_mll_met_2dpdf),
                                       ratio_realmll)

##################################################
# Load combined signal + background dataset
##################################################
combfile = ROOT.TFile.Open("signal_plus_background.root", "READ")
combtree = combfile.Get("event_tree")

weightVar = ROOT.RooRealVar("weight_nominal", "weight_nominal", -1, 1)
variables = ROOT.RooArgSet(mll, met, weightVar)
combdataset = ROOT.RooDataSet("combdataset", "Background plus signal MC",
                               variables,
                               RF.Import(combtree),
                               RF.WeightVar(weightVar))

print(f"\nDataset loaded: {combdataset.numEntries()} entries, "
      f"sum of weights = {combdataset.sumEntries():.2f}")

##################################################
# Extended 2D model: n_sig * sig_2dpdf + n_bkg * bkg_2dpdf
##################################################
n_total = combdataset.sumEntries()
n_sig = ROOT.RooRealVar("n_sig", "Signal yield",     n_total * 0.1, 0, n_total * 10) # made-up initial value
n_bkg = ROOT.RooRealVar("n_bkg", "Background yield", n_total * 0.9, 0, n_total * 10) # made-up initial value

model_2dpdf = ROOT.RooAddPdf("model_2dpdf", "Signal + Background 2D model",
                              ROOT.RooArgList(sigtot_mll_met_2dpdf, bkgtot_mll_met_2dpdf),
                              ROOT.RooArgList(n_sig, n_bkg))

##################################################
# Fit
##################################################
fit_result = model_2dpdf.fitTo(combdataset,
                                RF.Extended(True),
                                RF.Save(True),
                                RF.SumW2Error(True),
                                RF.PrintLevel(1))

##################################################
# Report
##################################################
print("\n=== Fit result ===")
fit_result.Print("v")
print(f"\nSignal yield:     {n_sig.getVal():.4f} +/- {n_sig.getError():.4f}")
print(f"Background yield: {n_bkg.getVal():.4f} +/- {n_bkg.getError():.4f}")
print(f"Fit status: {fit_result.status()}  (0 = OK)")
print(f"EDM:        {fit_result.edm():.3e}")

##################################################
# Save
##################################################
outfile = ROOT.TFile("fit2D_result_signal_plus_background.root", "RECREATE")
fit_result.Write("fit_result")
outfile.Close()
print("Saved to fit2D_result_signal_plus_background.root")

##################################################
# 1D projection plots (mll and MET)
# Each canvas: main pad (data + model components) + pull pad
##################################################
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
CMS.setCMSStyle()
CMS.SetExtraText("Simulation Preliminary")
CMS.cms_lumi = "13 TeV"
CMS.cms_energy = ""

obs_configs = [
    (mll, 30,   60,  120, "m(ll) [GeV]", "mll_bkgAndSig_fit_to_SplusB"),
    (met, 40,    0, 1200, "MET [GeV]",   "met_bkgAndSig_fit_to_SplusB"),
]

n_sig_val  = n_sig.getVal()
n_bkg_val  = n_bkg.getVal()
n_tot_val  = n_sig_val + n_bkg_val

for obs, nBins, xmin, xmax, xlabel, plotname in obs_configs:
    frame = obs.frame(RF.Bins(nBins), RF.Range(xmin, xmax), RF.Title(""))

    # Total model first (needed before pullHist can reference it)
    model_2dpdf.plotOn(frame,
                       RF.Name("total"),
                       RF.Normalization(n_tot_val, ROOT.RooAbsReal.NumEvent),
                       RF.LineColor( ROOT.TColor.GetColor("#9c9ca1") ), # grey
                       RF.LineWidth(2))

    # Signal component
    model_2dpdf.plotOn(frame,
                       RF.Components("sigtot_mll_met_2dpdf"),
                       RF.Name("signal"),
                       RF.Normalization(n_tot_val, ROOT.RooAbsReal.NumEvent),
                       RF.LineColor(ROOT.TColor.GetColor("#bd1f01")), # red
                       RF.LineStyle(ROOT.kDashed),
                       RF.LineWidth(2))

    # Background component
    model_2dpdf.plotOn(frame,
                       RF.Components("bkgtot_mll_met_2dpdf"),
                       RF.Name("background"),
                       RF.Normalization(n_tot_val, ROOT.RooAbsReal.NumEvent),
                       RF.LineColor(ROOT.TColor.GetColor("#5790fc")), # blue
                       RF.LineStyle(ROOT.kDashed),
                       RF.LineWidth(2))

    # Data on top (plotted last so markers appear above curves)
    combdataset.plotOn(frame,
                       RF.Binning(nBins),
                       RF.Name("data"),
                       RF.MarkerColor(ROOT.kBlack),
                       RF.LineColor(ROOT.kBlack),
                       RF.MarkerStyle(ROOT.kFullCircle),
                       RF.MarkerSize(0.8))

    # Pull histogram: (data - model) / sigma_data, bin-by-bin
    pull_hist = frame.pullHist("data", "total")
    pull_hist.SetMarkerStyle(ROOT.kFullCircle)
    pull_hist.SetMarkerSize(0.8)

    frame_pull = obs.frame(RF.Bins(nBins), RF.Range(xmin, xmax), RF.Title(""))
    frame_pull.addPlotable(pull_hist, "P")

    # Canvas split 70/30 between main and pull
    c = ROOT.TCanvas(plotname, plotname, 800, 800)
    pad_main = ROOT.TPad("pad_main", "", 0, 0.28, 1, 1.0)
    pad_pull = ROOT.TPad("pad_pull", "", 0, 0.00, 1, 0.28)
    pad_main.SetBottomMargin(0.02)
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
    frame.SetMaximum(1.6 * frame.GetMaximum())
    frame.Draw()

    leg = ROOT.TLegend(0.40, 0.65, 0.90, 0.90)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.042)
    leg.AddEntry(frame.findObject("data"),       "MC (background only)", "PE")
    leg.AddEntry(frame.findObject("total"),       "Total model", "L")
    leg.AddEntry(frame.findObject("signal"),      f"Signal (n_{{sig}} = {n_sig.getVal():.2f} +/- {n_sig.getError():.2f})", "L")
    leg.AddEntry(frame.findObject("background"),  f"Background (n_{{bkg}} = {n_bkg.getVal():.2f} +/- {n_bkg.getError():.2f})", "L")
    leg.Draw()
    # CMS_lumi draws the CMS text and luminosity information on the *specified* pad: https://cmsstyle.readthedocs.io/en/latest/reference/#cmsstyle.cmsstyle.CMS_lumi
    CMS.CMS_lumi(pad_main , iPosX=0)


    # info = ROOT.TLatex()
    # info.SetNDC()
    # info.SetTextSize(0.038)
    # info.DrawLatex(0.14, 0.88, f"Fit status: {fit_result.status()}   EDM: {fit_result.edm():.2e}")

    pad_pull.cd()
    frame_pull.GetYaxis().SetTitle("Pull")
    frame_pull.GetYaxis().CenterTitle(True)
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
    zero_line.SetLineColor(ROOT.kBlack)
    zero_line.SetLineWidth(1)
    zero_line.Draw("SAME")

    c.SaveAs(f"{plotname}.pdf")
    c.SaveAs(f"{plotname}.png")
    print(f"Created {plotname}.pdf / .png, copying to /eos/user/s/skkwan/www/higgsino/studies/mll-MET-fit-2D/overall-fit")
    os.system("mv *.pdf /eos/user/s/skkwan/www/higgsino/studies/mll-MET-fit-2D/overall-fit")
    os.system("mv *.png /eos/user/s/skkwan/www/higgsino/studies/mll-MET-fit-2D/overall-fit")

    # Log-scale version
    logname = f"{plotname}_logscale"
    c_log = ROOT.TCanvas(logname, logname, 800, 800)
    pad_main_log = ROOT.TPad("pad_main_log", "", 0, 0.28, 1, 1.0)
    pad_pull_log = ROOT.TPad("pad_pull_log", "", 0, 0.00, 1, 0.28)
    pad_main_log.SetBottomMargin(0.02)
    pad_main_log.SetLeftMargin(0.12)
    pad_pull_log.SetTopMargin(0.04)
    pad_pull_log.SetBottomMargin(0.35)
    pad_pull_log.SetLeftMargin(0.12)
    pad_main_log.Draw()
    pad_pull_log.Draw()

    pad_main_log.cd()
    pad_main_log.SetLogy()
    frame.GetYaxis().SetTitle("Events / bin")
    frame.SetMaximum(frame.GetMaximum() * 1e6)
    frame.SetMinimum(1e-10)
    frame.Draw()
    leg.Draw()

    # CMS_lumi draws the CMS text and luminosity information on the specified pad: https://cmsstyle.readthedocs.io/en/latest/reference/#cmsstyle.cmsstyle.CMS_lumi
    CMS.CMS_lumi(pad_main_log, iPosX=0)


    pad_pull_log.cd()
    frame_pull.Draw()
    zero_line.Draw("SAME")

    c_log.SaveAs(f"{logname}.pdf")
    c_log.SaveAs(f"{logname}.png")
    print(f"Created {logname}.pdf / .png, copying to /eos/user/s/skkwan/www/higgsino/studies/mll-MET-fit-2D/overall-fit")
    os.system("mv *.pdf /eos/user/s/skkwan/www/higgsino/studies/mll-MET-fit-2D/overall-fit")
    os.system("mv *.png /eos/user/s/skkwan/www/higgsino/studies/mll-MET-fit-2D/overall-fit")

sig_file = ROOT.TFile.Open("../pdf_fit/snapshot_TChiZH_650_1_SR_mll_MET_fit_scheme.root", "READ")
bkg_file = ROOT.TFile.Open("../pdf_fit/backgrounds_for_2D_fit.root", "READ")

def sum_weights(tree):
    return sum(getattr(ev, "weight_nominal") for ev in tree)

sig_tree = sig_file.Get("event_tree")
bkg_tree = bkg_file.Get("event_tree")

n_sig    = sig_tree.GetEntries()
n_bkg    = bkg_tree.GetEntries()
yield_sig = sum_weights(sig_tree)
yield_bkg = sum_weights(bkg_tree)

print(f"TChiZH signal events:    {n_sig}")
print(f"TChiZH signal yield:     {yield_sig:.4f}")
print(f"Background events:       {n_bkg}")
print(f"Background yield:        {yield_bkg:.4f}")
print(f"Total events:            {n_sig + n_bkg}")
print(f"Total yield:             {yield_sig + yield_bkg:.4f}")

sig_file.Close()
bkg_file.Close()