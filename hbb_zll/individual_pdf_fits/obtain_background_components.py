# Run in ROOT 6.38 (do not do cmsenv)

import os
import ROOT
from ROOT import RooFit as RF
import numpy as np
from array import array
import cmsstyle as CMS

def plotFit(name, rooVar, dataset, pdf, dataLabel, fitLabel, plotname,
            nFloatParams=2, outdir="", getOverflow=True, doLog=False, varBinEdges=None):
    """
    varBinEdges: optional list of bin edges used for chi^2 computation and data
                 display on the frame. The ratio panel uses the same variable bins
                 so numerator and denominator are on identical scales.
    """
    nBinsDisplay = 50
    xmin = rooVar.getMin()
    xmax = rooVar.getMax()

    # Build chi^2 binning
    if varBinEdges is not None:
        nBinsChi2 = len(varBinEdges) - 1
        chi2_binning = ROOT.RooBinning(nBinsChi2, array('d', varBinEdges))
        data_binning_arg = ROOT.RooFit.Binning(chi2_binning)
    else:
        data_binning_arg = ROOT.RooFit.Binning(nBinsDisplay)

    # Frame uses many bins so the PDF curve looks smooth
    frame = rooVar.frame(nBinsDisplay)

    # Plot data and PDF on the frame first so frame.GetMaximum() reflects the
    # actual displayed scale (PDF is normalized to data by RooFit, so
    # pdf.createHistogram() alone underestimates the true curve maximum).
    dataset.plotOn(frame, ROOT.RooFit.Name("data_" + plotname),
                   ROOT.RooFit.LineColor(ROOT.TColor.GetColor("#9c9ca1")),
                   ROOT.RooFit.LineWidth(2),
                   ROOT.RooFit.MarkerColor(ROOT.TColor.GetColor("#9c9ca1")),
                   ROOT.RooFit.MarkerSize(1),
                   data_binning_arg)

    pdf.plotOn(frame, ROOT.RooFit.Name("pdf_" + plotname),
               ROOT.RooFit.LineColor(ROOT.TColor.GetColor("#bd1f01")),
               ROOT.RooFit.LineWidth(2),
               ROOT.RooFit.LineStyle(1),
               ROOT.RooFit.MarkerSize(0),
               ROOT.RooFit.Binning(nBinsDisplay))

    roo_curve = frame.getCurve("pdf_" + plotname)
    chi2_per_ndf = frame.chiSquare("pdf_" + plotname, "data_" + plotname, nFloatParams)

    # Rescale to events/GeV: divide data by actual bin width and PDF by frame
    # reference bin width so both are on the same density axis.
    # Must happen after chi2 (which reads unscaled values) and before y_max.
    w_ref = (xmax - xmin) / nBinsDisplay
    roo_hist = frame.getHist("data_" + plotname)
    xs, ys = roo_hist.GetX(), roo_hist.GetY()
    ey_lo, ey_hi = roo_hist.GetEYlow(), roo_hist.GetEYhigh()
    for i in range(roo_hist.GetN()):
        if varBinEdges is not None:
            w_i = next(
                varBinEdges[j+1] - varBinEdges[j]
                for j in range(len(varBinEdges) - 1)
                if varBinEdges[j] <= xs[i] <= varBinEdges[j+1]
            )
        else:
            w_i = w_ref
        roo_hist.SetPoint(i, xs[i], ys[i] / w_i)
        roo_hist.SetPointEYlow(i, ey_lo[i] / w_i)
        roo_hist.SetPointEYhigh(i, ey_hi[i] / w_i)
    for i in range(roo_curve.GetN()):
        roo_curve.SetPoint(i, roo_curve.GetX()[i], roo_curve.GetY()[i] / w_ref)
    y_max = 1.5 * max(roo_hist.GetY()[i] + roo_hist.GetEYhigh()[i]
                      for i in range(roo_hist.GetN()))

    y_min = 0
    if doLog:
        y_min = 1e-10
        y_max = y_max * 1000

    # If you want to have "CMS (Private work)" and not the "CMS (Preliminary)" default text,
    # you need both of the following lines
    CMS.SetExtraText("Private work")
    CMS.SetCmsText("CMS", font=62, size=0.76)
    CMS.SetLumi(250, unit="fb", run="2018")

    canv = CMS.cmsDiCanvas("canv_" + plotname, x_min=xmin, x_max=xmax, y_min=y_min, y_max=y_max,
                           r_min=0, r_max=2,
                           nameXaxis=f"{name} / GeV",
                           nameYaxis="Events / GeV",
                           nameRatio="MC/Pred",
                           square=True, iPos=0)
    canv.SetRightMargin(0.05)
    CMS.UpdatePad(canv)

    canv.cd(1)
    if doLog:
        ROOT.gPad.SetLogy()
    CMS.UpdatePad(canv)

    leg = CMS.cmsLeg(0.2, 0.89 - 0.05 * 4, 0.95, 0.89, textSize=0.04)
    CMS.SetLumi(250, unit="fb", run="2018")

    leg.AddEntry(frame.findObject("data_" + plotname), dataLabel)
    leg.AddEntry(frame.findObject("pdf_" + plotname), fitLabel)
    leg.SetHeader(f"2018 SR: background MET ( #chi^{{2}}/ndf = {chi2_per_ndf:.2f})")
    frame.Draw("SAME")

    # Ratio panel
    canv.cd(2)
    if varBinEdges is not None:
        data_ratio = dataset.createHistogram("histo_ratio_" + plotname, rooVar,
                                             ROOT.RooFit.Binning(chi2_binning))
        pdf_ratio = pdf.createHistogram("hpdf_ratio_" + plotname, rooVar,
                                        ROOT.RooFit.Binning(chi2_binning))
        if pdf_ratio.Integral() > 0:
            pdf_ratio.Scale(data_ratio.Integral() / pdf_ratio.Integral())
        data_ratio.Divide(pdf_ratio)
    else:
        data_hist = dataset.createHistogram("histo_" + plotname, rooVar,
                                            ROOT.RooFit.Binning(nBinsDisplay, xmin, xmax))
        if getOverflow:
            nb = data_hist.GetNbinsX()
            data_hist.SetBinContent(nb, data_hist.GetBinContent(nb) + data_hist.GetBinContent(nb + 1))
        data_ratio = data_hist.Clone()
        prediction = data_hist.Clone()
        for i in range(1, data_ratio.GetNbinsX() + 1):
            prediction.SetBinContent(i, roo_curve.Eval(data_ratio.GetBinCenter(i)) * w_ref)
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

    # Workaround for "CMS Preliminary" being constantly cut off
    canv.SaveAs(f"{fname}.eps")
    os.system(f"gs -q -dBATCH -dNOPAUSE -dSAFER -dEPSCrop -dPDFSETTINGS=/prepress -sDEVICE=pdfwrite -dEmbedAllFonts=true -dSubsetFonts=true -sOutputFile={fname}.pdf {fname}.eps && rm {fname}.eps")

    ROOT.gStyle.SetImageScaling(3.)
    canv.SaveAs(f"{fname}.png")
    ROOT.gStyle.SetImageScaling(1.)
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

    # Plot all objects on the frame first so frame.GetMaximum() reflects the
    # actual displayed scale after RooFit normalization.
    peaking_dataset.plotOn(frame, ROOT.RooFit.Name("peaking_data_" + plotname),
                           ROOT.RooFit.LineColor(ROOT.TColor.GetColor("#5790fc")),
                           ROOT.RooFit.LineWidth(2),
                           ROOT.RooFit.MarkerColor(ROOT.TColor.GetColor("#5790fc")),
                           ROOT.RooFit.MarkerSize(0.5),
                           ROOT.RooFit.Binning(nBins))

    nonpeak_dataset.plotOn(frame, ROOT.RooFit.Name("nonpeak_data_" + plotname),
                           ROOT.RooFit.LineColor(ROOT.TColor.GetColor("#964a8b")),
                           ROOT.RooFit.LineWidth(2),
                           ROOT.RooFit.MarkerColor(ROOT.TColor.GetColor("#964a8b")),
                           ROOT.RooFit.MarkerSize(0.5),
                           ROOT.RooFit.Binning(nBins))

    peaking_pdf.plotOn(frame, ROOT.RooFit.Name("peaking_pdf_" + plotname),
                       ROOT.RooFit.LineColor(ROOT.TColor.GetColor("#f89c20")),
                       ROOT.RooFit.LineWidth(2),
                       ROOT.RooFit.LineStyle(1),
                       ROOT.RooFit.MarkerSize(0),
                       ROOT.RooFit.Normalization(peaking_dataset.sumEntries(), ROOT.RooAbsReal.NumEvent),
                       ROOT.RooFit.Binning(nBins))

    nonpeak_pdf.plotOn(frame, ROOT.RooFit.Name("nonpeak_pdf_" + plotname),
                       ROOT.RooFit.LineColor(ROOT.TColor.GetColor("#e42536")),
                       ROOT.RooFit.LineWidth(2),
                       ROOT.RooFit.LineStyle(1),
                       ROOT.RooFit.MarkerSize(0),
                       ROOT.RooFit.Normalization(nonpeak_dataset.sumEntries(), ROOT.RooAbsReal.NumEvent),
                       ROOT.RooFit.Binning(nBins))

    y_min = 0
    y_max = 1.5 * frame.GetMaximum()
    if doLog:
        y_min = 1e-10
        y_max = y_max * 1000

    canv = CMS.cmsDiCanvas("canv_" + plotname, x_min=xmin, x_max=xmax, y_min=y_min, y_max=y_max,
                           r_min=0, r_max=2,
                           nameXaxis="MET / GeV",
                           nameYaxis="Shape (A.U.)",
                           nameRatio="MC/Pred",
                           square=True, iPos=0)
    canv.SetRightMargin(0.05)
    CMS.UpdatePad(canv)

    pad1 = canv.cd(1)
    pad1.SetTopMargin(0.15)
    if doLog:
        ROOT.gPad.SetLogy()
    CMS.UpdatePad(canv)

    leg = CMS.cmsLeg(0.2, 0.89 - 0.05 * 4, 0.9, 0.89, textSize=0.04)
    leg.SetHeader("2018 SR: background MET components")
    CMS.SetLumi(250, unit="fb", run="2018")

    leg.AddEntry(frame.findObject("peaking_data_" + plotname), "Peaking-in-m(ll) background")
    leg.AddEntry(frame.findObject("nonpeak_data_" + plotname), "Non-peaking-in-m(ll) background")
    leg.AddEntry(frame.findObject("peaking_pdf_" + plotname), "Peaking Gumbel fit")
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


def plotInputHistogram(name, rooVar, dataset, label, plotname, color,
                       varBinEdges=None, nBins=50, outdir="", doLog=False):
    """
    Plot the raw input histogram with variable or uniform binning. From Claude
    """
    xmin = rooVar.getMin()
    xmax = rooVar.getMax()

    if varBinEdges is not None:
        binning = ROOT.RooBinning(len(varBinEdges) - 1, array('d', varBinEdges))
        h = dataset.createHistogram("h_input_" + plotname, rooVar,
                                    ROOT.RooFit.Binning(binning))
    else:
        h = dataset.createHistogram("h_input_" + plotname, rooVar,
                                    ROOT.RooFit.Binning(nBins, xmin, xmax))

    n_events = dataset.sumEntries()
    h.SetLineColor(ROOT.TColor.GetColor(color))
    h.SetMarkerColor(ROOT.TColor.GetColor(color))
    h.SetMarkerStyle(ROOT.kFullCircle)
    h.SetMarkerSize(0.2)
    h.SetLineWidth(2)

    y_min = 0
    y_max = h.GetMaximum() * 1.5
    if doLog:
        y_min = 1e-10
        y_max = h.GetMaximum() * 1000

    CMS.SetCmsText("CMS", font=62, size=0.76)
    CMS.SetLumi(250, run="2018")
    canv = CMS.cmsCanvas("canv_input_" + plotname, xmin, xmax, y_min, y_max,
                         nameXaxis=f"{name} / GeV",
                         nameYaxis="Events",
                         square=True, iPos=0) # extraSpace=0.05, yTitOffset=1.6, iPos=0)
    canv.SetRightMargin(0.1)
    CMS.UpdatePad(canv)

    if doLog:
        ROOT.gPad.SetLogy()
        CMS.UpdatePad(canv)

    leg = CMS.cmsLeg(0.25, 0.92 - 0.05 * 3, 0.90, 0.92, textSize=0.03)
    leg.SetMargin(0.12)
    # leg.SetBorderSize(1)
    # leg.SetLineColor(ROOT.kBlack)
    # leg.SetLineWidth(1)
    leg.AddEntry(h, label + f" (n = {n_events:.2f})")

    CMS.cmsObjectDraw(h, "E1")
    CMS.cmsObjectDraw(leg)
 
    if varBinEdges is not None:
        for obj in canv.GetListOfPrimitives():
            if obj.InheritsFrom("TH1"):
                obj.GetXaxis().SetNdivisions(5, False)
                break
        ROOT.gPad.Modified()
        ROOT.gPad.Update()

    CMS.UpdatePad(canv)

    fname = "input_" + plotname + ("-log" if doLog else "")
    canv.SaveAs(f"{fname}.eps")
    os.system(f"gs -q -dBATCH -dNOPAUSE -dSAFER -dEPSCrop -dPDFSETTINGS=/prepress -sDEVICE=pdfwrite -dEmbedAllFonts=true -dSubsetFonts=true -sOutputFile={fname}.pdf {fname}.eps && rm {fname}.eps")

    ROOT.gStyle.SetImageScaling(3.)
    canv.SaveAs(f"{fname}.png")
    ROOT.gStyle.SetImageScaling(1.)
    if outdir:
        os.system(f"mv {fname}.* {outdir}")

    del canv
    del leg


def plotMETPDFsOnly(rooVar, peaking_pdf, nonpeak_pdf, plotname, outdir="", getOverflow=True, doLog=False):
    nBins = 50
    xmin = rooVar.getMin()
    xmax = rooVar.getMax()

    frame = rooVar.frame(nBins)

    # Plot PDFs on the frame first so frame.GetMaximum() reflects the actual
    # displayed scale after RooFit normalization.
    peaking_pdf.plotOn(frame, ROOT.RooFit.Name("peaking_pdf_" + plotname),
                       ROOT.RooFit.LineColor(ROOT.TColor.GetColor("#f89c20")),
                       ROOT.RooFit.LineWidth(2),
                       ROOT.RooFit.LineStyle(1),
                       ROOT.RooFit.MarkerSize(0),
                       ROOT.RooFit.Normalization(1.0, ROOT.RooAbsReal.NumEvent),
                       ROOT.RooFit.Binning(nBins))

    nonpeak_pdf.plotOn(frame, ROOT.RooFit.Name("nonpeak_pdf_" + plotname),
                       ROOT.RooFit.LineColor(ROOT.TColor.GetColor("#e42536")),
                       ROOT.RooFit.LineWidth(2),
                       ROOT.RooFit.LineStyle(1),
                       ROOT.RooFit.MarkerSize(0),
                       ROOT.RooFit.Normalization(1.0, ROOT.RooAbsReal.NumEvent),
                       ROOT.RooFit.Binning(nBins))

    y_min = 0
    y_max = 1.5 * frame.GetMaximum()
    if doLog:
        y_min = 1e-10
        y_max = y_max * 1000

    canv = CMS.cmsDiCanvas("canv_" + plotname, x_min=xmin, x_max=xmax, y_min=y_min, y_max=y_max,
                           r_min=0, r_max=2,
                           nameXaxis="MET / GeV",
                           nameYaxis="Shape (A.U.)",
                           nameRatio="MC/Pred",
                           square=True, iPos=0)
    canv.SetRightMargin(0.05)
    CMS.UpdatePad(canv)

    canv.cd(1)
    if doLog:
        ROOT.gPad.SetLogy()
    CMS.UpdatePad(canv)

    leg = CMS.cmsLeg(0.2, 0.89 - 0.05 * 2, 0.9, 0.89, textSize=0.04)
    leg.SetHeader("2018 SR: background MET components")
    CMS.SetLumi(250, unit="fb", run="2018")
    leg.AddEntry(frame.findObject("peaking_pdf_" + plotname), "Peaking Gumbel fit")
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
##### Variable bin boundaries for MET chi^2
##################################################
met_var_bins = [200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 500, 800, 1200]

ROOT.gStyle.SetImageScaling(1.)  # set to 3 only immediately before PNG saves

##################################################
##### Define fit observables
##################################################
mll = ROOT.RooRealVar("m_ll", "m_ll", 60, 120)
met = ROOT.RooRealVar("met", "met", 200, 1200)

##################################################
###### Retrieve background dataset (prepared with reformat.py)
##################################################
bkg_peaking_filepath = '../input_data/backgrounds_peaking.root'
bkg_nonpeak_filepath = '../input_data/backgrounds_nonpeak.root'
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


os.makedirs("individual_fit_results", exist_ok=True)
f = ROOT.TFile("individual_fit_results/fitresult_background_all_except_ZPeak.root", "RECREATE")
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
               "bkg_nonpeak_met_gumbel", varBinEdges=met_var_bins, doLog=doLog)
    plotFit("MET", met, bkg_peaking_dataset_met, bkgpeaking_met,
               "Peaking-in-m(ll) background",
               f"Gumbel fit (#mu={mu_peaking_met.getVal():.1f} #pm {mu_peaking_met.getError():.1f}, b={b_peaking_met.getVal():.1f}#pm{b_peaking_met.getError():.1f})",
               "bkg_peaking_met_gumbel", varBinEdges=met_var_bins, doLog=doLog)
    plotFit("m(ll)", mll, bkg_nonpeak_dataset_mll, bkgnonpeak_mll,
               "Non-peaking-in-m(ll) background",
               f"Exponential fit (a={a_nonpeak_mll.getVal():.2f}  #pm {a_nonpeak_mll.getError():.2f})",
               "bkg_nonpeak_mll_exponential", doLog=doLog)
    plotInputHistogram("MET", met, bkg_nonpeak_dataset_met,
                       "Non-peaking-in-m(ll) background", "bkg_nonpeak_met_input",
                       color="#964a8b", varBinEdges=met_var_bins, doLog=doLog)
    plotInputHistogram("MET", met, bkg_peaking_dataset_met,
                       "Peaking-in-m(ll) background", "bkg_peaking_met_input",
                       color="#5790fc", varBinEdges=met_var_bins, doLog=doLog)
    plotInputHistogram("m(ll)", mll, bkg_nonpeak_dataset_mll,
                       "Non-peaking-in-m(ll) background", "bkg_nonpeak_mll_input",
                       color="#964a8b", nBins=50, doLog=doLog)
# plotMETPDFTogether(met,
#                    bkg_peaking_dataset, bkg_nonpeak_dataset,
#                    bkg_peaking_met_pdf, bkgnonpeak_met,
#                    "bkg_met_components_gumbel")
# plotMETPDFsOnly(met,
#                 bkg_peaking_met_pdf, bkgnonpeak_met,
#                 "bkg_met_pdfs_only_gumbel")



os.system("mv *.png /eos/user/s/skkwan/www/higgsino/studies/mll-MET-fit-2D/background_shapes")
os.system("mv *.pdf /eos/user/s/skkwan/www/higgsino/studies/mll-MET-fit-2D/background_shapes")

