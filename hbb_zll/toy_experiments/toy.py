import ROOT
import cmsstyle as CMS
import os
from uncertainties import ufloat

def get_signal_model(m1=650, m2=1):
    """
    Get the signal model for a given mass point.
    """
    # Load results
    signalresultsfile = ROOT.TFile.Open(f"../individual_pdf_fits/individual_fit_results/fitresult_signal_{m1}_{m2}.root", "READ")

    workspace              = signalresultsfile.Get(f"workspace_{m1}_{m2}")
    sig_result             = signalresultsfile.Get("sig_result")

    # sigtot_mll_met_2dpdf = sig_dcb_mll (mll) x pdf_of_spline (MET)
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

    components = [
        signalresultsfile, 
        workspace,
        sig_result, mean_mll, sigmal_mll, sigmar_mll, alphal_mll, alphar_mll, nl_mll, nr_mll, 
        sig_dcb_mll, pdf_of_spline, 
    ]

    return sigtot_mll_met_2dpdf, components


def get_background_model():
    """
    Get the total background model, leaving the ratio r (peaking/total) floating.
    """
    # Load results
    bkgresultsfile    = ROOT.TFile.Open("../individual_pdf_fits/individual_fit_results/fitresult_background_all_except_ZPeak.root", "READ")
    zpeakresultsfile  = ROOT.TFile.Open("../zpeak_fit/initial_zPeak_fit_result.root", "READ")

    bkg_nonpeak_result_met = bkgresultsfile.Get("bkg_nonpeak_result_met")
    bkg_nonpeak_result_mll = bkgresultsfile.Get("bkg_nonpeak_result_mll")
    bkg_peaking_result_met = bkgresultsfile.Get("bkg_peaking_result_met")
    zpeak_result           = zpeakresultsfile.Get("zPeak_CRZ_fit_result")

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

    # Leave the ratio_peaking floating
    ratio_peaking = ROOT.RooRealVar("ratio_peaking", "ratio_peaking", 0.088, 0, 1)

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

    bkgtot_mll_met_2dpdf = ROOT.RooAddPdf("bkgtot_mll_met_2dpdf", "bkgtot_mll_met_2dpdf",
                                           ROOT.RooArgList(bkgpeaking_mll_met_2dpdf,
                                                           bkgnonpeak_mll_met_2dpdf),
                                           ratio_peaking)

    components = [
        bkgresultsfile, zpeakresultsfile,
        bkg_nonpeak_result_met, bkg_nonpeak_result_mll,
        bkg_peaking_result_met, zpeak_result,
        bkgnonpeak_met, bkgnonpeak_mll, bkgnonpeak_mll_met_2dpdf,
        bkgpeaking_met, bkgpeaking_mll, bkgpeaking_mll_met_2dpdf,
        ratio_peaking,
    ]
    return bkgtot_mll_met_2dpdf, components

def make_toy_plot(frame, obs, nBins, xmin, xmax, xlabel, plotname, doLog,
                  n_sig, n_bkg, n_peak_val, n_peak_err, n_nonpeak_val, n_nonpeak_err,
                  ratio_val, mass_point, eos_dir="."):
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

    frame_pull = obs.frame(ROOT.RooFit.Bins(nBins), ROOT.RooFit.Range(xmin, xmax), ROOT.RooFit.Title(""))
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
    leg.AddEntry(frame.findObject("data"),        "Toy data", "PE")
    leg.AddEntry(frame.findObject("total"),       "Total model", "L")
    # leg.AddEntry(frame.findObject("signal"),      f"Signal ({m1}, {m2}) GeV (n_{{sig}} = {n_sig.getVal():.2f} +/- {n_sig.getError():.2f})", "L")
    leg.AddEntry(frame.findObject("bkg_peaking"), f"Peaking background (n_{{peak}} = {n_peak_val:.2f} +/- {n_peak_err:.2f})", "L")
    leg.AddEntry(frame.findObject("bkg_nonpeak"), f"Non-peaking background (n_{{nonpeak}} = {n_nonpeak_val:.2f} +/- {n_nonpeak_err:.2f})", "L")
    r_dummy = ROOT.TLine()
    r_dummy.SetLineWidth(0)
    r_dummy.SetLineColor(0)
    leg.AddEntry(r_dummy, f"n_{{bkg}} = {n_bkg.getVal():.2f} +/- {n_bkg.getError():.2f}", "L")
    leg.AddEntry(r_dummy, f"Ratio r (peaking / total) = {ratio_val:.3f}", "L")

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

    canv.SaveAs(f"{plotname}.pdf")
    canv.SaveAs(f"{plotname}.png")
    print(f"Created {plotname}.pdf / .png")
    if eos_dir != ".":
        os.system(f"mv {plotname}.pdf {eos_dir}/")
        os.system(f"mv {plotname}.png {eos_dir}/")


def plot_toy_fit(toys, model, result, met, mll, mass_point,
                 plotname_prefix="toy_fit", eos_dir="."):
    """
    Plot 1D projections of the toy dataset with the best-fit model overlaid,
    in the same style as fit2D.py.
    """
    n_sig = result.floatParsFinal().find("n_sig")
    n_bkg = result.floatParsFinal().find("n_bkg")
    ratio_peaking = result.floatParsFinal().find("ratio_peaking")

    # n_sig_val     = n_sig.getVal()
    n_bkg_val     = n_bkg.getVal()
    n_bkg_err     = n_bkg.getError()
    # n_tot_val     = n_sig_val + n_bkg_val
    n_tot_val = n_bkg_val # TODO: fix this
    n_peak_val    = n_bkg_val * ratio_peaking.getVal()
    n_nonpeak_val = n_bkg_val * (1 - ratio_peaking.getVal())
    n_peak_err    = n_bkg_err * ratio_peaking.getVal()
    n_nonpeak_err = n_bkg_err * (1 - ratio_peaking.getVal())

    obs_configs = [
        (mll, 40,  60.,  120., "m(ll) [GeV]", f"mll_{plotname_prefix}"),
        (met, 60, 200., 1200., "MET [GeV]",   f"met_{plotname_prefix}"),
    ]

    for obs, nBins, xmin, xmax, xlabel, plotname in obs_configs:
        frame = obs.frame(ROOT.RooFit.Bins(nBins), ROOT.RooFit.Range(xmin, xmax), ROOT.RooFit.Title(""))

        model.plotOn(frame,
                     ROOT.RooFit.Name("total"),
                     ROOT.RooFit.Normalization(n_tot_val, ROOT.RooAbsReal.NumEvent),
                     ROOT.RooFit.LineColor(ROOT.TColor.GetColor("#9c9ca1")),
                     ROOT.RooFit.LineWidth(2))

        # model.plotOn(frame,
        #              ROOT.RooFit.Components("sigtot_mll_met_2dpdf"),
        #              ROOT.RooFit.Name("signal"),
        #              ROOT.RooFit.Normalization(n_tot_val, ROOT.RooAbsReal.NumEvent),
        #              ROOT.RooFit.LineColor(ROOT.TColor.GetColor("#bd1f01")),
        #              ROOT.RooFit.LineStyle(ROOT.kDashed),
        #              ROOT.RooFit.LineWidth(2))

        model.plotOn(frame,
                     ROOT.RooFit.Components("bkgpeaking_mll_met_2dpdf"),
                     ROOT.RooFit.Name("bkg_peaking"),
                     ROOT.RooFit.Normalization(n_tot_val, ROOT.RooAbsReal.NumEvent),
                     ROOT.RooFit.LineColor(ROOT.TColor.GetColor("#964a8b")),
                     ROOT.RooFit.LineStyle(ROOT.kDashed),
                     ROOT.RooFit.LineWidth(2))

        model.plotOn(frame,
                     ROOT.RooFit.Components("bkgnonpeak_mll_met_2dpdf"),
                     ROOT.RooFit.Name("bkg_nonpeak"),
                     ROOT.RooFit.Normalization(n_tot_val, ROOT.RooAbsReal.NumEvent),
                     ROOT.RooFit.LineColor(ROOT.TColor.GetColor("#3f90da")),
                     ROOT.RooFit.LineStyle(ROOT.kDashed),
                     ROOT.RooFit.LineWidth(2))

        toys.plotOn(frame,
                    ROOT.RooFit.Binning(nBins),
                    ROOT.RooFit.Name("data"),
                    ROOT.RooFit.MarkerColor(ROOT.kBlack),
                    ROOT.RooFit.LineColor(ROOT.kBlack),
                    ROOT.RooFit.MarkerStyle(ROOT.kFullCircle),
                    ROOT.RooFit.MarkerSize(0.8))

        plot_kwargs = dict(obs=obs, nBins=nBins, xmin=xmin, xmax=xmax, xlabel=xlabel,
                           plotname=plotname,
                           n_sig=n_sig, n_bkg=n_bkg,
                           n_peak_val=n_peak_val, n_peak_err=n_peak_err,
                           n_nonpeak_val=n_nonpeak_val, n_nonpeak_err=n_nonpeak_err,
                           ratio_val=ratio_peaking.getVal(), mass_point=mass_point,
                           eos_dir="/eos/user/s/skkwan/www/higgsino/studies/mll-MET-fit-2D/toys")

        make_toy_plot(frame, doLog=False, **plot_kwargs)
        make_toy_plot(frame, doLog=True,  **plot_kwargs)


if __name__ == "__main__" :
    m1 = 650
    m2 = 1

    # Observables
    met = ROOT.RooRealVar("met", "met", 200, 1200)
    mll = ROOT.RooRealVar("m_ll", "m_ll", 60, 120)

    # Generate background and signal toys in the same function
    n_toys = 200

    # Get the background plus signal model
    n_sig_in = 0
    n_bkg_in = 100
    n_sig = ROOT.RooRealVar("n_sig", "n_sig", n_sig_in, -1, n_toys*5)
    n_bkg = ROOT.RooRealVar("n_bkg", "n_bkg", n_bkg_in, 0, n_toys*5)

    bkg_model, bkg_components = get_background_model()
    sig_model, sig_components = get_signal_model(m1, m2)
    # model = ROOT.RooAddPdf("total_pdf", "total_pdf",
    #                              ROOT.RooArgList(sig_model, bkg_model),
    #                              ROOT.RooArgList(n_sig, n_bkg))
    model = ROOT.RooAddPdf("total_pdf", "total_pdf",
                                 ROOT.RooArgList(bkg_model),
                                 ROOT.RooArgList(n_bkg))                             
    ROOT.RooRandom.randomGenerator().SetSeed(0)
    toys = model.generate(ROOT.RooArgSet(met, mll), n_toys)

    result = model.fitTo(toys, ROOT.RooFit.Save(True))

    # Plot 
    met_frame = met.frame()
    toys.plotOn(met_frame)
    model.plotOn(met_frame)
    model.plotOn(met_frame, ROOT.RooFit.Components("total_pdf"), ROOT.RooFit.LineStyle(ROOT.kDashed))
    canv = ROOT.TCanvas("canv", "canv", 800, 600)
    met_frame.Draw()
    canv.SaveAs("test.pdf")
    canv.SaveAs("test.png")
    os.system("mv test.* /eos/user/s/skkwan/www/higgsino/studies/mll-MET-fit-2D/toys")


    n_sig.Print()
    n_bkg.Print()

    # print("\n=== Fit result ===")
    # result.Print("v")
    # print(f"Fit status: {result.status()}  (0 = OK)")
    # print(f"EDM:        {result.edm():.3e}")

    # nb_comb_fit = ufloat(result.floatParsFinal().find("n_bkg").getVal(), result.floatParsFinal().find("n_bkg").getError())
    # ns_comb_fit = ufloat(result.floatParsFinal().find("n_sig").getVal(), result.floatParsFinal().find("n_sig").getError())

    # print(f"{nb_comb_fit:.2f}")
    # print(f"{ns_comb_fit:.2f}")

    plot_toy_fit(toys, model, result, met, mll,
                 mass_point=(m1, m2))

    # #  Pull:
    # n_sig_pull = (result.floatParsFinal().find("n_sig").getVal() - n_sig_in) / result.floatParsFinal().find("n_sig").getError()
    # print(f"Signal pull: {n_sig_pull}")