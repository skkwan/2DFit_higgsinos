# Run in ROOT 6.38 (do not do cmsenv)

import os
import ROOT
from ROOT import RooFit as RF
import numpy as np
from array import array
import cmsstyle as CMS

# From Claude and modified
def plotFit(name, rooVar, dataset, pdf, dataLabel, fitLabel, plotname, nFloatParams=2, outdir="", getOverflow=True, doLog=False):
    nBins = 50
    xmin = rooVar.getMin()
    xmax = rooVar.getMax()

    frame = rooVar.frame(nBins)

    leg = CMS.cmsLeg(0.3, 0.89 - 0.05 * 4, 0.95, 0.89, textSize=0.04)
    CMS.SetLumi("")

    data_hist = dataset.createHistogram("histo_" + plotname, rooVar, ROOT.RooFit.Binning(nBins, xmin, xmax))
    if getOverflow:
        data_hist.SetBinContent(data_hist.GetNbinsX(), data_hist.GetBinContent(data_hist.GetNbinsX()) + data_hist.GetBinContent(data_hist.GetNbinsX() + 1))
    pdf_hist = pdf.createHistogram("hpdf_" + plotname, rooVar, ROOT.RooFit.Binning(nBins))
    if getOverflow:
        pdf_hist.SetBinContent(pdf_hist.GetNbinsX(), pdf_hist.GetBinContent(pdf_hist.GetNbinsX()) + pdf_hist.GetBinContent(pdf_hist.GetNbinsX() + 1))
    y_min = 0
    y_max = 1.8 * max(data_hist.GetMaximum(), pdf_hist.GetMaximum())
    if doLog:
        y_min = 1e-10
        y_max = y_max * 1000

    canv = CMS.cmsDiCanvas("canv_" + plotname, x_min=xmin, x_max=xmax, y_min=y_min, y_max=y_max,
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
    chi2_per_ndf = frame.chiSquare("pdf_" + plotname, "data_" + plotname, nFloatParams)
    leg.SetHeader(f"2018 SR: background MET ( #chi^{{2}}/ndf = {chi2_per_ndf:.2f})")
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


def plotMETPDFTogether(rooVar, peaking_dataset, nonpeak_dataset,
                       peaking_pdf, nonpeak_pdf, plotname, outdir="", getOverflow=True, doLog=False):
    nBins = 50
    xmin = rooVar.getMin()
    xmax = rooVar.getMax()

    frame = rooVar.frame(nBins)

    leg = CMS.cmsLeg(0.3, 0.89 - 0.05 * 4, 0.9, 0.89, textSize=0.04)
    leg.SetHeader("2018 SR: background MET components")
    CMS.SetLumi("")

    peak_hist = peaking_dataset.createHistogram("histo_peak_" + plotname, rooVar, ROOT.RooFit.Binning(nBins, xmin, xmax))
    nonpeak_hist = nonpeak_dataset.createHistogram("histo_nonpeak_" + plotname, rooVar, ROOT.RooFit.Binning(nBins, xmin, xmax))
    if getOverflow:
        peak_hist.SetBinContent(peak_hist.GetNbinsX(), peak_hist.GetBinContent(peak_hist.GetNbinsX()) + peak_hist.GetBinContent(peak_hist.GetNbinsX() + 1))
        nonpeak_hist.SetBinContent(nonpeak_hist.GetNbinsX(), nonpeak_hist.GetBinContent(nonpeak_hist.GetNbinsX()) + nonpeak_hist.GetBinContent(nonpeak_hist.GetNbinsX() + 1))
    peak_pdf_hist = peaking_pdf.createHistogram("hpdf_peak_" + plotname, rooVar, ROOT.RooFit.Binning(nBins))
    nonpeak_pdf_hist = nonpeak_pdf.createHistogram("hpdf_nonpeak_" + plotname, rooVar, ROOT.RooFit.Binning(nBins))
    if getOverflow:
        peak_pdf_hist.SetBinContent(peak_pdf_hist.GetNbinsX(), peak_pdf_hist.GetBinContent(peak_pdf_hist.GetNbinsX()) + peak_pdf_hist.GetBinContent(peak_pdf_hist.GetNbinsX() + 1))
        nonpeak_pdf_hist.SetBinContent(nonpeak_pdf_hist.GetNbinsX(), nonpeak_pdf_hist.GetBinContent(nonpeak_pdf_hist.GetNbinsX()) + nonpeak_pdf_hist.GetBinContent(nonpeak_pdf_hist.GetNbinsX() + 1))
    y_min = 0
    y_max = 1.8 * max(peak_hist.GetMaximum(), nonpeak_hist.GetMaximum(),
                      peak_pdf_hist.GetMaximum(), nonpeak_pdf_hist.GetMaximum())
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

    peaking_dataset.plotOn(frame, ROOT.RooFit.Name("peaking_data_" + plotname),
                           ROOT.RooFit.LineColor(ROOT.TColor.GetColor("#5790fc")),
                           ROOT.RooFit.LineWidth(2),
                           ROOT.RooFit.MarkerColor(ROOT.TColor.GetColor("#5790fc")),
                           ROOT.RooFit.MarkerSize(1),
                           ROOT.RooFit.Binning(nBins))
    leg.AddEntry(frame.findObject("peaking_data_" + plotname), "Peaking-in-m(ll) background")

    nonpeak_dataset.plotOn(frame, ROOT.RooFit.Name("nonpeak_data_" + plotname),
                           ROOT.RooFit.LineColor(ROOT.TColor.GetColor("#964a8b")),
                           ROOT.RooFit.LineWidth(2),
                           ROOT.RooFit.MarkerColor(ROOT.TColor.GetColor("#964a8b")),
                           ROOT.RooFit.MarkerSize(1),
                           ROOT.RooFit.Binning(nBins))
    leg.AddEntry(frame.findObject("nonpeak_data_" + plotname), "Non-peaking-in-m(ll) background")

    peaking_pdf.plotOn(frame, ROOT.RooFit.Name("peaking_pdf_" + plotname),
                       ROOT.RooFit.LineColor(ROOT.TColor.GetColor("#f89c20")),
                       ROOT.RooFit.LineWidth(2),
                       ROOT.RooFit.LineStyle(1),
                       ROOT.RooFit.MarkerSize(0),
                       ROOT.RooFit.Normalization(peaking_dataset.sumEntries(), ROOT.RooAbsReal.NumEvent),
                       ROOT.RooFit.Binning(nBins))
    leg.AddEntry(frame.findObject("peaking_pdf_" + plotname), "Peaking Gumbel fit")

    nonpeak_pdf.plotOn(frame, ROOT.RooFit.Name("nonpeak_pdf_" + plotname),
                       ROOT.RooFit.LineColor(ROOT.TColor.GetColor("#e42536")),
                       ROOT.RooFit.LineWidth(2),
                       ROOT.RooFit.LineStyle(1),
                       ROOT.RooFit.MarkerSize(0),
                       ROOT.RooFit.Normalization(nonpeak_dataset.sumEntries(), ROOT.RooAbsReal.NumEvent),
                       ROOT.RooFit.Binning(nBins))
    leg.AddEntry(frame.findObject("nonpeak_pdf_" + plotname), "Non-peaking Gumbel fit")

    frame.Draw("SAME")

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


def plotMETPDFsOnly(rooVar, peaking_pdf, nonpeak_pdf, plotname, outdir="", getOverflow=True, doLog=False):
    nBins = 50
    xmin = rooVar.getMin()
    xmax = rooVar.getMax()

    frame = rooVar.frame(nBins)

    leg = CMS.cmsLeg(0.3, 0.89 - 0.05 * 2, 0.9, 0.89, textSize=0.04)
    leg.SetHeader("2018 SR: background MET components")
    CMS.SetLumi("")

    peak_pdf_hist = peaking_pdf.createHistogram("hpdf_peak_" + plotname, rooVar, ROOT.RooFit.Binning(nBins))
    nonpeak_pdf_hist = nonpeak_pdf.createHistogram("hpdf_nonpeak_" + plotname, rooVar, ROOT.RooFit.Binning(nBins))
    if getOverflow:
        peak_pdf_hist.SetBinContent(peak_pdf_hist.GetNbinsX(), peak_pdf_hist.GetBinContent(peak_pdf_hist.GetNbinsX()) + peak_pdf_hist.GetBinContent(peak_pdf_hist.GetNbinsX() + 1))
        nonpeak_pdf_hist.SetBinContent(nonpeak_pdf_hist.GetNbinsX(), nonpeak_pdf_hist.GetBinContent(nonpeak_pdf_hist.GetNbinsX()) + nonpeak_pdf_hist.GetBinContent(nonpeak_pdf_hist.GetNbinsX() + 1))
    y_min = 0
    y_max = 1.8 * max(peak_pdf_hist.GetMaximum(), nonpeak_pdf_hist.GetMaximum())
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

    peaking_pdf.plotOn(frame, ROOT.RooFit.Name("peaking_pdf_" + plotname),
                       ROOT.RooFit.LineColor(ROOT.TColor.GetColor("#f89c20")),
                       ROOT.RooFit.LineWidth(2),
                       ROOT.RooFit.LineStyle(1),
                       ROOT.RooFit.MarkerSize(0),
                       ROOT.RooFit.Normalization(1.0, ROOT.RooAbsReal.NumEvent),
                       ROOT.RooFit.Binning(nBins))
    leg.AddEntry(frame.findObject("peaking_pdf_" + plotname), "Peaking Gumbel fit")

    nonpeak_pdf.plotOn(frame, ROOT.RooFit.Name("nonpeak_pdf_" + plotname),
                       ROOT.RooFit.LineColor(ROOT.TColor.GetColor("#e42536")),
                       ROOT.RooFit.LineWidth(2),
                       ROOT.RooFit.LineStyle(1),
                       ROOT.RooFit.MarkerSize(0),
                       ROOT.RooFit.Normalization(1.0, ROOT.RooAbsReal.NumEvent),
                       ROOT.RooFit.Binning(nBins))
    leg.AddEntry(frame.findObject("nonpeak_pdf_" + plotname), "Non-peaking Gumbel fit")

    frame.Draw("SAME")

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
# Non-peaking: MET and m(ll)
###########################################################################
# 1(a) First do NON-peaking background by itself: MET
bkg_nonpeak_dataset_met = ROOT.RooDataSet("bkg_nonpeak_dataset_met", "bkg_nonpeak_dataset_met", ROOT.RooArgSet(met, weightXyear), ROOT.RooFit.Import(bkg_nonpeak_tree), ROOT.RooFit.WeightVar(weightXyear))
mu_nonpeak_met = ROOT.RooRealVar('mu_nonpeak_met', 'mu_nonpeak_met', 225, 100, 300)
b_nonpeak_met = ROOT.RooRealVar('b_nonpeak_met', 'b_nonpeak_met', 60, 20, 100)
bkgnonpeak_met = ROOT.RooGenericPdf("bkgnonpeak_met", "bkgnonpeak_met", "1/b_nonpeak_met * exp(-(@0 - mu_nonpeak_met)/b_nonpeak_met - exp(-(@0 - mu_nonpeak_met)/b_nonpeak_met))",
                        ROOT.RooArgList(met, mu_nonpeak_met, b_nonpeak_met))
# Fit:
bkg_nonpeak_result_met = bkgnonpeak_met.fitTo(bkg_nonpeak_dataset_met, RF.Save(), RF.SumW2Error(True))
bkg_nonpeak_params_met = bkgnonpeak_met.getParameters(bkg_nonpeak_dataset_met)
print(bkg_nonpeak_params_met.Print("v"))

# 1(b) NON-peaking background by itself: m(ll)
bkg_nonpeak_dataset_mll = ROOT.RooDataSet("bkg_nonpeak_dataset_mll", "bkg_nonpeak_dataset_mll", ROOT.RooArgSet(mll, weightXyear), ROOT.RooFit.Import(bkg_nonpeak_tree), ROOT.RooFit.WeightVar(weightXyear))
a_nonpeak_mll = ROOT.RooRealVar("a_nonpeak_mll", "a_nonpeak_mll", -0.03, -1, 1)
bkgnonpeak_mll = ROOT.RooExponential("bkgnonpeak_mll", "bkgnonpeak_mll", mll, a_nonpeak_mll)
# Fit 
bkg_nonpeak_result_mll = bkgnonpeak_mll.fitTo(bkg_nonpeak_dataset_mll, RF.Save(), RF.SumW2Error(True))
bkg_nonpeak_params_mll = bkgnonpeak_mll.getParameters(bkg_nonpeak_dataset_mll)

###########################################################################
# Peaking: MET (m(ll) already gotten from a different fit)
###########################################################################
# 2(a) PEAKING background by itself: MET 
mu_peaking_met = ROOT.RooRealVar('mu_peaking_met', 'mu_peaking_met', 50, 30, 400)   # adequate value somewhere around 246 based on hand-drawn plots
b_peaking_met = ROOT.RooRealVar('b_peaking_met', 'b_peaking_met', 40, 20, 100)    # adequate value somewhere around 41 based on hand-drawn plots
bkgpeaking_met = ROOT.RooGenericPdf("bkgpeaking_met", "bkgpeaking_met", "1/b_peaking_met * exp(-(@0 - mu_peaking_met)/b_peaking_met - exp(-(@0 - mu_peaking_met)/b_peaking_met))",
                        ROOT.RooArgList(met, mu_peaking_met, b_peaking_met))
bkg_peaking_dataset_met = ROOT.RooDataSet("bkg_peaking_dataset_met", "bkg_peaking_dataset_met", ROOT.RooArgSet(met, weightXyear), ROOT.RooFit.Import(bkg_peaking_tree), ROOT.RooFit.WeightVar(weightXyear))
# Fit 
bkg_peaking_result_met = bkgpeaking_met.fitTo(bkg_peaking_dataset_met, RF.Save(), RF.SumW2Error(True))
bkg_peaking_params_met = bkgpeaking_met.getParameters(bkg_peaking_dataset_met)


f = ROOT.TFile("fitresult_background_all_except_ZPeak.root", "RECREATE")
bkg_nonpeak_result_met.Write("bkg_nonpeak_result_met")
bkg_nonpeak_result_mll.Write("bkg_nonpeak_result_mll")
bkg_peaking_result_met.Write("bkg_peaking_result_met")

bkg_nonpeak_params_met.Write("bkg_nonpeak_params_met")
bkg_nonpeak_params_mll.Write("bkg_nonpeak_params_mll")
bkg_peaking_params_met.Write("bkg_peaking_params_met")

f.Close()


for doLog in [True, False]:
    plotFit("MET", met, bkg_nonpeak_dataset_met, bkgnonpeak_met,
               "Non-peaking-in-m(ll) background",
               f"Gumbel fit (#mu={mu_nonpeak_met.getVal():.1f} #pm {mu_nonpeak_met.getError():.1f}, b={b_nonpeak_met.getVal():.1f}#pm{b_nonpeak_met.getError():.1f})",
               "bkg_nonpeak_met_gumbel", doLog=doLog)
    plotFit("MET", met, bkg_peaking_dataset_met, bkgpeaking_met,
               "Peaking-in-m(ll) background",
               f"Gumbel fit (#mu={mu_peaking_met.getVal():.1f} #pm {mu_peaking_met.getError():.1f}, b={b_peaking_met.getVal():.1f}#pm{b_peaking_met.getError():.1f})",
               "bkg_peaking_met_gumbel", doLog=doLog)
    plotFit("m(ll)", mll, bkg_nonpeak_dataset_mll, bkgnonpeak_mll,
               "Non-peaking-in-m(ll) background",
               f"Exponential fit (a={a_nonpeak_mll.getVal():.2f}  #pm {a_nonpeak_mll.getError():.2f})",
               "bkg_nonpeak_mll_exponential", doLog=doLog)
# plotMETPDFTogether(met,
#                    bkg_peaking_dataset, bkg_nonpeak_dataset,
#                    bkg_peaking_met_pdf, bkgnonpeak_met,
#                    "bkg_met_components_gumbel")
# plotMETPDFsOnly(met,
#                 bkg_peaking_met_pdf, bkgnonpeak_met,
#                 "bkg_met_pdfs_only_gumbel")



os.system("mv *.png /eos/user/s/skkwan/www/higgsino/studies/mll-MET-fit-2D/background_shapes")
os.system("mv *.pdf /eos/user/s/skkwan/www/higgsino/studies/mll-MET-fit-2D/background_shapes")
