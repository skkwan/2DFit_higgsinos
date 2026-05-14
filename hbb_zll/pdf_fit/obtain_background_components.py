# Run in ROOT 6.38 (do not do cmsenv)

import os
import ROOT
from ROOT import RooFit as RF
import numpy as np
from array import array
import cmsstyle as CMS

doLog = True

# From Claude
def plotMETFit(rooVar, dataset, pdf, dataLabel, fitLabel, plotname, outdir=""):
    nBins = 50
    xmin = rooVar.getMin()
    xmax = rooVar.getMax()

    frame = rooVar.frame(nBins)

    leg = CMS.cmsLeg(0.3, 0.89 - 0.05 * 4, 0.9, 0.89, textSize=0.04)
    leg.SetHeader("2018 signal region MC: background MET")
    CMS.SetLumi("")

    data_hist = dataset.createHistogram("histo_" + plotname, rooVar, ROOT.RooFit.Binning(nBins, xmin, xmax))
    y_min = 0
    y_max = 1.8 * data_hist.GetMaximum()
    if doLog:
        y_min = 1e-10
        y_max = y_max * 1000

    canv = CMS.cmsDiCanvas("canv_" + plotname, x_min=xmin, x_max=xmax, y_min=y_min, y_max=y_max,
                           r_min=0, r_max=2,
                           nameXaxis="MET / GeV",
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
                   ROOT.RooFit.LineColor(ROOT.TColor.GetColor("#9c9ca1")),
                   ROOT.RooFit.LineWidth(2),
                   ROOT.RooFit.MarkerColor(ROOT.TColor.GetColor("#9c9ca1")),
                   ROOT.RooFit.MarkerSize(1),
                   ROOT.RooFit.Binning(nBins))
    leg.AddEntry(frame.findObject("data_" + plotname), dataLabel)

    pdf.plotOn(frame, ROOT.RooFit.Name("pdf_" + plotname),
               ROOT.RooFit.LineColor(ROOT.TColor.GetColor("#bd1f01")),
               ROOT.RooFit.LineWidth(2),
               ROOT.RooFit.LineStyle(1),
               ROOT.RooFit.MarkerSize(0),
               ROOT.RooFit.Binning(nBins))
    leg.AddEntry(frame.findObject("pdf_" + plotname), fitLabel)

    roo_curve = frame.getCurve("pdf_" + plotname)
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
##### Define fit observables 
##################################################
mll = ROOT.RooRealVar("m_ll", "m_ll", 60, 120)
met = ROOT.RooRealVar("met", "met", 200, 1200)

##################################################
###### Retrieve background dataset (prepared with reformat.py)
##################################################
bkg_peaking_filepath = 'backgrounds_peaking.root'
bkg_nonpeak_filepath = 'backgrounds_nonpeak.root'
bkg_peaking_file = ROOT.TFile.Open(bkg_peaking_filepath, "READ")
bkg_nonpeak_file = ROOT.TFile.Open(bkg_nonpeak_filepath, "READ")

bkg_peaking_tree = bkg_peaking_file.Get("event_tree")
bkg_nonpeak_tree = bkg_nonpeak_file.Get("event_tree")

weightXyear = ROOT.RooRealVar("weight_nominal", "weight_nominal", -1, 1)
variables = ROOT.RooArgSet(mll, met, weightXyear)

###########################################################################
# MET: Declare the background 1D PDFs and multiply them into 2D PDFs
###########################################################################
# First do NON-peaking background by itself
bkg_nonpeak_dataset = ROOT.RooDataSet("bkg_nonpeak_dataset", "bkg_nonpeak_dataset", ROOT.RooArgSet(met, weightXyear), ROOT.RooFit.Import(bkg_nonpeak_tree), ROOT.RooFit.WeightVar(weightXyear))
# Background fake mll component in MET dimension
a_nonpeak_met = ROOT.RooRealVar('a_nonpeak_met', 'a_nonpeak_met', -0.005, -0.1, -0.0001)
bkg_nonpeak_met_pdf = ROOT.RooExponential("bkgfakemll_met", "bkgfakemll_met", met, a_nonpeak_met)
# Fit and get results
bkg_nonpeak_result = bkg_nonpeak_met_pdf.fitTo(bkg_nonpeak_dataset, RF.Save(), RF.SumW2Error(True))
bkg_nonpeak_params = bkg_nonpeak_met_pdf.getParameters(bkg_nonpeak_dataset)
print(bkg_nonpeak_params.Print("v"))

# Then do PEAKING background by itself
bkg_peaking_dataset = ROOT.RooDataSet("bkg_peaking_dataset", "bkg_peaking_dataset", ROOT.RooArgSet(met, weightXyear), ROOT.RooFit.Import(bkg_peaking_tree), ROOT.RooFit.WeightVar(weightXyear))
a_peaking_met = ROOT.RooRealVar('a_peaking_met', 'a_peaking_met', -0.005, -0.1, -0.0001)
bkg_peaking_met_pdf = ROOT.RooExponential("bkgrealmll_met", "bkgrealmll_met", met, a_peaking_met)
# Fit
bkg_peaking_result = bkg_peaking_met_pdf.fitTo(bkg_peaking_dataset, RF.Save(), RF.SumW2Error(True))
bkg_peaking_params = bkg_peaking_met_pdf.getParameters(bkg_peaking_dataset)
print(bkg_peaking_params.Print("v"))

# Then do TOTAL background
bkg_total_file = ROOT.TFile.Open('backgrounds_for_2D_fit.root', "READ")
bkg_total_tree = bkg_total_file.Get("event_tree")
bkg_total_dataset = ROOT.RooDataSet("bkg_total_dataset", "bkg_total_dataset", ROOT.RooArgSet(met, weightXyear), ROOT.RooFit.Import(bkg_total_tree), ROOT.RooFit.WeightVar(weightXyear))
a_total_met = ROOT.RooRealVar('a_total_met', 'a_total_met', -0.005, -0.1, -0.0001)
bkg_total_met_pdf = ROOT.RooExponential("bkg_total_met", "bkg_total_met", met, a_total_met)
# Fit
bkg_total_result = bkg_total_met_pdf.fitTo(bkg_total_dataset, RF.Save(), RF.SumW2Error(True))
bkg_total_params = bkg_total_met_pdf.getParameters(bkg_total_dataset)
print(bkg_total_params.Print("v"))

f = ROOT.TFile("fitresult_background_MET.root", "RECREATE")
bkg_nonpeak_result.Write("bkg_nonpeak_result")
bkg_peaking_result.Write("bkg_peaking_result")
bkg_total_result.Write("bkg_total_result")
f.Close()

plotMETFit(met, bkg_nonpeak_dataset, bkg_nonpeak_met_pdf,
           "Non-peaking-in-m(ll) background",
           f"Exponential (a = {a_nonpeak_met.getVal():.4f} #pm {a_nonpeak_met.getError():.4f})",
           "bkg_nonpeak_met_exp")
plotMETFit(met, bkg_peaking_dataset, bkg_peaking_met_pdf,
           "Peaking-in-m(ll) background",
           f"Exponential (a = {a_peaking_met.getVal():.4f} #pm {a_peaking_met.getError():.4f})",
           "bkg_peaking_met_exp")
plotMETFit(met, bkg_total_dataset, bkg_total_met_pdf,
           "Total background",
           f"Exponential (a = {a_total_met.getVal():.4f} #pm {a_total_met.getError():.4f})",
           "bkg_total_met_exp")

os.system("mv *.png /eos/user/s/skkwan/www/higgsino/studies/mll-MET-fit-2D/background_shapes")
os.system("mv *.pdf /eos/user/s/skkwan/www/higgsino/studies/mll-MET-fit-2D/background_shapes")

