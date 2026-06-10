# fit2D.py
# Run in ROOT 6.38 (do not do cmsenv)
#
# Performs a 2D (mll, MET) extended maximum-likelihood fit using a
# signal + background model for a given signal mass point.
#
# hasSignal=False: fit to background-only MC (backgrounds_total.root)
# hasSignal=True:  fit to signal + background MC (signal_plus_backgrounds_<m1>_<m2>.root)
#
# All shape parameters are fixed from pre-fit results.
# The only floating parameters are n_sig, n_bkg, and ratio_peaking.

import ROOT
from ROOT import RooFit as RF
import cmsstyle as CMS
import os
import numpy as np


def make_plot(frame, obs, nBins, xmin, xmax, xlabel, plotname, doLog,
              n_sig, n_bkg, n_peak_val, n_peak_err, n_nonpeak_val, n_nonpeak_err,
              ratio_val, ratio_err, mass_point, data_legend_label, eos_dir):
    if doLog:
        plotname = f"{plotname}_logscale"
        y_min = 1e-6
        y_max = frame.GetMaximum() * 1e10
        frame.SetMinimum(y_min)
        frame.SetMaximum(y_max)
    else:
        y_min = 0
        y_max = 2.2 * frame.GetMaximum()

    pull_hist = frame.pullHist("data", "total")
    pull_hist.SetMarkerStyle(ROOT.kFullCircle)
    pull_hist.SetMarkerSize(0.8)

    frame_pull = obs.frame(RF.Bins(nBins), RF.Range(xmin, xmax), RF.Title(""))
    frame_pull.addPlotable(pull_hist, "P")

    CMS.SetExtraText("Private work")
    CMS.SetCmsText("CMS", font=62, size=0.76)
    CMS.SetLumi(250, unit="fb", run="2018")

    canv = CMS.cmsDiCanvas("canv_" + plotname, x_min=xmin, x_max=xmax, y_min=y_min, y_max=y_max,
                           r_min=-2, r_max=2,
                           nameXaxis=f"{xlabel}",
                           nameYaxis="Events",
                           nameRatio="Pull",
                           square=True, extraSpace=0.01, iPos=0.)
    canv.SetRightMargin(0.05)
    CMS.UpdatePad(canv)

    canv.cd(1)
    if doLog:
        ROOT.gPad.SetLogy()

    m1, m2 = mass_point
    leg = ROOT.TLegend(0.25, 0.50, 0.90, 0.90)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.035)
    leg.AddEntry(frame.findObject("data"),       data_legend_label, "PE")
    leg.AddEntry(frame.findObject("total"),      "Total model", "L")
    leg.AddEntry(frame.findObject("signal"),     f"Signal ({m1}, {m2}) GeV (n_{{sig}} = {n_sig.getVal():.2f} +/- {n_sig.getError():.2f})", "L")
    leg.AddEntry(frame.findObject("bkg_peaking"), f"Peaking background (n_{{peak}} = {n_peak_val:.2f} +/- {n_peak_err:.2f})", "L")
    leg.AddEntry(frame.findObject("bkg_nonpeak"), f"Non-peaking background (n_{{nonpeak}} = {n_nonpeak_val:.2f} +/- {n_nonpeak_err:.2f})", "L")
    r_dummy = ROOT.TLine()
    r_dummy.SetLineWidth(0)
    r_dummy.SetLineColor(0)
    leg.AddEntry(r_dummy, f"Total background yield: n_{{bkg}} = {n_bkg.getVal():.2f} +/- {n_bkg.getError():.2f}", "L")
    leg.AddEntry(r_dummy, f"Ratio r (peaking/ total) = {ratio_val:.2f}  #pm {ratio_err:.2f}", "L")

    frame.Draw("SAME")

    canv.cd(2)
    frame_pull.Draw("SAME")
    zero_line = ROOT.TLine(xmin, 0, xmax, 0)
    zero_line.SetLineColor(ROOT.kBlack)
    zero_line.SetLineWidth(1)
    zero_line.SetLineStyle(ROOT.kDashed)
    zero_line.Draw("SAME")

    canv.cd(1)
    CMS.cmsObjectDraw(leg)
    CMS.UpdatePad(canv)

    canv.SaveAs(f"{plotname}.eps")
    os.system(f"gs -q -dBATCH -dNOPAUSE -dSAFER -dEPSCrop -dPDFSETTINGS=/prepress -sDEVICE=pdfwrite "
              f"-dEmbedAllFonts=true -dSubsetFonts=true -sOutputFile={plotname}.pdf {plotname}.eps && rm {plotname}.eps")
    canv.SaveAs(f"{plotname}.png")
    print(f"Created {plotname}.pdf / .png, copying to {eos_dir}")
    os.system(f"mv *.pdf {eos_dir}")
    os.system(f"mv *.png {eos_dir}")


def run_fit(hasSignal, mass_point):
    m1, m2 = mass_point
    mp_str = f"{m1}_{m2}"

    ##################################################
    # Load pre-fit results
    ##################################################
    signalresultsfile = ROOT.TFile.Open(f"../individual_pdf_fits/individual_fit_results/fitresult_signal_{mp_str}.root", "READ")
    bkgresultsfile    = ROOT.TFile.Open("../individual_pdf_fits/individual_fit_results/fitresult_background_all_except_ZPeak.root", "READ")
    zpeakresultsfile  = ROOT.TFile.Open("../zpeak_fit/initial_zPeak_fit_result.root", "READ")

    workspace              = signalresultsfile.Get(f"workspace_{m1}_{m2}")
    sig_result             = signalresultsfile.Get("sig_result")
    bkg_nonpeak_result_met = bkgresultsfile.Get("bkg_nonpeak_result_met")
    bkg_nonpeak_result_mll = bkgresultsfile.Get("bkg_nonpeak_result_mll")
    bkg_peaking_result_met = bkgresultsfile.Get("bkg_peaking_result_met")
    zpeak_result           = zpeakresultsfile.Get("zPeak_CRZ_fit_result")

    ##################################################
    # Observables
    # Use the workspace's met variable so that pdf_of_spline
    # (which was built against it) shares the same object.
    ##################################################
    met = workspace.var("met")
    mll = ROOT.RooRealVar("m_ll", "m_ll", 60, 120)

    ##################################################
    # Signal 2D PDF (fixed shape)
    # sigtot_mll_met_2dpdf = sig_dcb_mll (mll) x pdf_of_spline (MET)
    ##################################################
    mean_mll   = sig_result.floatParsFinal().find(f"mean_mll_{m1}_{m2}")
    sigmal_mll = sig_result.floatParsFinal().find(f"sigmal_mll_{m1}_{m2}")
    sigmar_mll = sig_result.floatParsFinal().find(f"sigmar_mll_{m1}_{m2}")
    alphal_mll = sig_result.floatParsFinal().find(f"alphal_mll_{m1}_{m2}")
    alphar_mll = sig_result.floatParsFinal().find(f"alphar_mll_{m1}_{m2}")
    nl_mll     = sig_result.floatParsFinal().find(f"nl_mll_{m1}_{m2}")
    nr_mll     = sig_result.floatParsFinal().find(f"nr_mll_{m1}_{m2}")
    for v in [mean_mll, sigmal_mll, sigmar_mll, alphal_mll, alphar_mll, nl_mll, nr_mll]:
        v.setConstant(True)

    sig_dcb_mll = ROOT.RooCrystalBall("sig_dcb_mll", "sig_dcb_mll",
                                       mll, mean_mll, sigmal_mll, sigmar_mll,
                                       alphal_mll, nl_mll, alphar_mll, nr_mll)

    pdf_of_spline = workspace.pdf(f"pdf_of_spline_{m1}_{m2}")

    sigtot_mll_met_2dpdf = ROOT.RooProdPdf("sigtot_mll_met_2dpdf", "sigtot_mll_met_2dpdf",
                                            ROOT.RooArgList(sig_dcb_mll, pdf_of_spline))

    ##################################################
    # Background 2D PDF (fixed shape)
    # bkgtot_mll_met_2dpdf = ratio_peaking * bkgpeaking_2dpdf
    #                      + (1 - ratio_peaking) * bkgnonpeak_2dpdf
    ##################################################

    # Non-peaking in m(ll)
    mu_nonpeak_met = bkg_nonpeak_result_met.floatParsFinal().find("mu_nonpeak_met")
    b_nonpeak_met  = bkg_nonpeak_result_met.floatParsFinal().find("b_nonpeak_met")
    a_nonpeak_mll  = bkg_nonpeak_result_mll.floatParsFinal().find("a_nonpeak_mll")
    for v in [mu_nonpeak_met, b_nonpeak_met, a_nonpeak_mll]:
        v.setConstant(True)

    bkgnonpeak_met = ROOT.RooGenericPdf(
        "bkgnonpeak_met", "bkgnonpeak_met",
        "1/b_nonpeak_met * exp(-(@0 - mu_nonpeak_met)/b_nonpeak_met"
        " - exp(-(@0 - mu_nonpeak_met)/b_nonpeak_met))",
        ROOT.RooArgList(met, mu_nonpeak_met, b_nonpeak_met))

    bkgnonpeak_mll = ROOT.RooExponential("bkgnonpeak_mll", "bkgnonpeak_mll",
                                          mll, a_nonpeak_mll)

    bkgnonpeak_mll_met_2dpdf = ROOT.RooProdPdf("bkgnonpeak_mll_met_2dpdf",
                                                "bkgnonpeak_mll_met_2dpdf",
                                                ROOT.RooArgList(bkgnonpeak_mll, bkgnonpeak_met))

    # Peaking in m(ll)
    mu_peaking_met = bkg_peaking_result_met.floatParsFinal().find("mu_peaking_met")
    b_peaking_met  = bkg_peaking_result_met.floatParsFinal().find("b_peaking_met")
    for v in [mu_peaking_met, b_peaking_met]:
        v.setConstant(True)

    bkgpeaking_met = ROOT.RooGenericPdf(
        "bkgpeaking_met", "bkgpeaking_met",
        "1/b_peaking_met * exp(-(@0 - mu_peaking_met)/b_peaking_met"
        " - exp(-(@0 - mu_peaking_met)/b_peaking_met))",
        ROOT.RooArgList(met, mu_peaking_met, b_peaking_met))

    zpeak_mean_mll   = zpeak_result.floatParsFinal().find("peak_mean_mll")
    zpeak_sigmal_mll = zpeak_result.floatParsFinal().find("peak_sigmal_mll")
    zpeak_sigmar_mll = zpeak_result.floatParsFinal().find("peak_sigmar_mll")
    zpeak_alphal_mll = zpeak_result.floatParsFinal().find("peak_alphal_mll")
    zpeak_nl_mll     = zpeak_result.floatParsFinal().find("peak_nl_mll")
    zpeak_alphar_mll = zpeak_result.floatParsFinal().find("peak_alphar_mll")
    zpeak_nr_mll     = zpeak_result.floatParsFinal().find("peak_nr_mll")
    for v in [zpeak_mean_mll, zpeak_sigmal_mll, zpeak_sigmar_mll,
              zpeak_alphal_mll, zpeak_nl_mll, zpeak_alphar_mll, zpeak_nr_mll]:
        v.setConstant(True)

    bkgpeaking_mll = ROOT.RooCrystalBall("bkgpeaking_mll", "bkgpeaking_mll",
                                          mll, zpeak_mean_mll, zpeak_sigmal_mll, zpeak_sigmar_mll,
                                          zpeak_alphal_mll, zpeak_nl_mll, zpeak_alphar_mll, zpeak_nr_mll)

    bkgpeaking_mll_met_2dpdf = ROOT.RooProdPdf("bkgpeaking_mll_met_2dpdf",
                                                "bkgpeaking_mll_met_2dpdf",
                                                ROOT.RooArgList(bkgpeaking_mll, bkgpeaking_met))

    # Total background: ratio between peaking and non-peaking is floating
    ratio_peaking = ROOT.RooRealVar("ratio_peaking", "ratio_peaking", 0.1, 0, 1)

    bkgtot_mll_met_2dpdf = ROOT.RooAddPdf("bkgtot_mll_met_2dpdf", "bkgtot_mll_met_2dpdf",
                                           ROOT.RooArgList(bkgpeaking_mll_met_2dpdf,
                                                           bkgnonpeak_mll_met_2dpdf),
                                           ratio_peaking)

    ##################################################
    # Dataset and per-run configuration
    ##################################################
    if hasSignal:
        input_file        = f"../input_data/signal_plus_backgrounds_{mp_str}.root"
        dataset_label     = "Background plus signal MC"
        data_legend_label = "MC (signal + background)"
        plot_prefix       = f"bkgAndSig_fit_to_SplusB_{mp_str}"
        eos_dir           = "/eos/user/s/skkwan/www/higgsino/studies/mll-MET-fit-2D/two_dimensional_fit"
        outfile_name      = f"two_dimensional_fit_results/fit2D_result_bkgAndSig_{mp_str}.root"
    else:
        input_file        = "../input_data/backgrounds_total.root"
        dataset_label     = "Background MC only"
        data_legend_label = "MC (background only)"
        plot_prefix       = f"bkgOnly_fit_to_SplusB_{mp_str}"
        eos_dir           = "/eos/user/s/skkwan/www/higgsino/studies/mll-MET-fit-2D/two_dimensional_fit"
        outfile_name      = f"two_dimensional_fit_results/fit2D_result_bkgOnly_{mp_str}.root"

    combfile = ROOT.TFile.Open(input_file, "READ")
    combtree = combfile.Get("event_tree")

    weightVar = ROOT.RooRealVar("weight_nominal", "weight_nominal", -1, 1)
    variables = ROOT.RooArgSet(mll, met, weightVar)
    combdataset = ROOT.RooDataSet("combdataset", dataset_label,
                                   variables,
                                   RF.Import(combtree),
                                   RF.WeightVar(weightVar))

    print(f"\nDataset loaded: {combdataset.numEntries()} entries, "
          f"sum of weights = {combdataset.sumEntries():.2f}")

    ##################################################
    # Extended 2D model: n_sig * sig_2dpdf + n_bkg * bkg_2dpdf
    ##################################################
    n_total = combdataset.sumEntries()
    n_sig = ROOT.RooRealVar("n_sig", "Signal yield",     n_total * 0.1, 0, n_total * 10)
    n_bkg = ROOT.RooRealVar("n_bkg", "Background yield", n_total * 0.9, 0, n_total * 10)

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
    os.makedirs("two_dimensional_fit_results", exist_ok=True)
    outfile = ROOT.TFile(outfile_name, "RECREATE")
    fit_result.Write("fit_result")
    outfile.Close()
    print(f"Saved to {outfile_name}")

    ##################################################
    # 1D projection plots (mll and MET)
    ##################################################
    n_sig_val     = n_sig.getVal()
    n_bkg_val     = n_bkg.getVal()
    n_bkg_err     = n_bkg.getError()
    n_tot_val     = n_sig_val + n_bkg_val
    ratio_val     = ratio_peaking.getVal()
    ratio_err     = ratio_peaking.getError()
    n_peak_val    = n_bkg_val * ratio_val
    n_nonpeak_val = n_bkg_val * (1 - ratio_val)
    n_peak_err    = np.sqrt((ratio_val * n_bkg_err)**2 + (n_bkg_val * ratio_err)**2)
    n_nonpeak_err = np.sqrt(((1 - ratio_val) * n_bkg_err)**2 + (n_bkg_val * ratio_err)**2)

    obs_configs = [
        (mll, 40,  60.,  120., "m(ll) [GeV]", f"mll_{plot_prefix}"),
        (met, 60,   0., 1200., "MET [GeV]",   f"met_{plot_prefix}"),
    ]

    for obs, nBins, xmin, xmax, xlabel, plotname in obs_configs:
        frame = obs.frame(RF.Bins(nBins), RF.Range(xmin, xmax), RF.Title(""))

        model_2dpdf.plotOn(frame,
                           RF.Name("total"),
                           RF.Normalization(n_tot_val, ROOT.RooAbsReal.NumEvent),
                           RF.LineColor(ROOT.TColor.GetColor("#9c9ca1")),
                           RF.LineWidth(2))

        model_2dpdf.plotOn(frame,
                           RF.Components("sigtot_mll_met_2dpdf"),
                           RF.Name("signal"),
                           RF.Normalization(n_tot_val, ROOT.RooAbsReal.NumEvent),
                           RF.LineColor(ROOT.TColor.GetColor("#bd1f01")),
                           RF.LineStyle(ROOT.kDashed),
                           RF.LineWidth(2))

        model_2dpdf.plotOn(frame,
                           RF.Components("bkgpeaking_mll_met_2dpdf"),
                           RF.Name("bkg_peaking"),
                           RF.Normalization(n_tot_val, ROOT.RooAbsReal.NumEvent),
                           RF.LineColor(ROOT.TColor.GetColor("#964a8b")),
                           RF.LineStyle(ROOT.kDashed),
                           RF.LineWidth(2))

        model_2dpdf.plotOn(frame,
                           RF.Components("bkgnonpeak_mll_met_2dpdf"),
                           RF.Name("bkg_nonpeak"),
                           RF.Normalization(n_tot_val, ROOT.RooAbsReal.NumEvent),
                           RF.LineColor(ROOT.TColor.GetColor("#3f90da")),
                           RF.LineStyle(ROOT.kDashed),
                           RF.LineWidth(2))

        combdataset.plotOn(frame,
                           RF.Binning(nBins),
                           RF.Name("data"),
                           RF.MarkerColor(ROOT.kBlack),
                           RF.LineColor(ROOT.kBlack),
                           RF.MarkerStyle(ROOT.kFullCircle),
                           RF.MarkerSize(0.8))

        plot_kwargs = dict(obs=obs, nBins=nBins, xmin=xmin, xmax=xmax, xlabel=xlabel,
                           plotname=plotname, n_sig=n_sig, n_bkg=n_bkg,
                           n_peak_val=n_peak_val, n_peak_err=n_peak_err,
                           n_nonpeak_val=n_nonpeak_val, n_nonpeak_err=n_nonpeak_err,
                           ratio_val=ratio_val, ratio_err=ratio_err,
                           mass_point=mass_point,
                           data_legend_label=data_legend_label,
                           eos_dir=eos_dir)

        make_plot(frame, doLog=False, **plot_kwargs)
        make_plot(frame, doLog=True,  **plot_kwargs)

    ##################################################
    # Event / yield summary (signal+background datasets only)
    ##################################################
    if hasSignal:
        def sum_weights(tree):
            return sum(getattr(ev, "weight_nominal") for ev in tree)

        sig_file = ROOT.TFile.Open(f"../input_data/snapshot_TChiZH_{mp_str}_SR_mll_MET_fit_scheme.root", "READ")
        bkg_file = ROOT.TFile.Open("../input_data/backgrounds_total.root", "READ")
        sig_tree = sig_file.Get("event_tree")
        bkg_tree = bkg_file.Get("event_tree")

        n_sig_entries = sig_tree.GetEntries()
        n_bkg_entries = bkg_tree.GetEntries()
        yield_sig     = sum_weights(sig_tree)
        yield_bkg     = sum_weights(bkg_tree)

        print(f"TChiZH ({m1}, {m2}) signal events:  {n_sig_entries}")
        print(f"TChiZH ({m1}, {m2}) signal yield:   {yield_sig:.4f}")
        print(f"Background events:                   {n_bkg_entries}")
        print(f"Background yield:                    {yield_bkg:.4f}")
        print(f"Total events:                        {n_sig_entries + n_bkg_entries}")
        print(f"Total yield:                         {yield_sig + yield_bkg:.4f}")

        sig_file.Close()
        bkg_file.Close()


if __name__ == "__main__":
    mass_points = [(650, 1), (500, 375)]
    for mp in mass_points:
        print(f"\n{'='*60}")
        print(f"  Background-only fit  |  mass point {mp}")
        print(f"{'='*60}")
        run_fit(hasSignal=False, mass_point=mp)

        print(f"\n{'='*60}")
        print(f"  Signal + background fit  |  mass point {mp}")
        print(f"{'='*60}")
        run_fit(hasSignal=True, mass_point=mp)
