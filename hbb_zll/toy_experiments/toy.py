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


def get_background_model(r=0.088):
    """
    Get the total background model, using a ratio r for the ratio of the peaking background to the total background.
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

    # Total background: ratio between peaking and non-peaking is set to r
    ratio_peaking = ROOT.RooRealVar("ratio_peaking", "ratio_peaking", 0.1, 0, 1)
    ratio_peaking.setVal(r)
    ratio_peaking.setConstant(True)

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

def gen_toys(met, mll, num_toys, model):
    """
    Generate (met, mll) toys for a given model and the given number of toys
    """
    toys = model.generate(ROOT.RooArgSet(met, mll), num_toys)
    return toys

def combine_signal_plus_background_toys(bkg_toys, sig_toys):
    """
    Combine signal and background toys and return the sum as its own object
    """
    total_toys = bkg_toys.Clone()
    total_toys.append(sig_toys) 
    return total_toys

def get_signal_plus_background_model(n_sig, n_bkg, m1=650, m2=1, r=0.088):
    """
    Get the signal plus background model for a specific signal mass point and ratio r (peaking background / total background).
    Uses n_sig and n_bkg as the signal and background yields.
    """
    bkg_model, bkg_components = get_background_model(r)
    sig_model, sig_components = get_signal_model(m1, m2)
    n_sig_var = ROOT.RooRealVar("n_sig", "n_sig", n_sig, 0, n_sig * 20)
    n_bkg_var = ROOT.RooRealVar("n_bkg", "n_bkg", n_bkg, 0, n_bkg * 20)
    total_model = ROOT.RooAddPdf("sig_bkg_2dpdf", "sig_bkg_2dpdf",
                                 ROOT.RooArgList(bkg_model, sig_model),
                                 ROOT.RooArgList(n_bkg_var, n_sig_var))
    total_components = bkg_components + sig_components + [bkg_model, sig_model, n_sig_var, n_bkg_var]
    return total_model, total_components

def do_combined_fit(model, toys):
    """
    Fit toys to model and analyze results.
    """
    result = model.fitTo(toys, ROOT.RooFit.Save(), ROOT.RooFit.SumW2Error(False))
    print("\n=== Fit result ===")
    result.Print("v")
    print(f"Fit status: {result.status()}  (0 = OK)")
    print(f"EDM:        {result.edm():.3e}")
    return result

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
    leg.AddEntry(frame.findObject("signal"),      f"Signal ({m1}, {m2}) GeV (n_{{sig}} = {n_sig.getVal():.2f} +/- {n_sig.getError():.2f})", "L")
    leg.AddEntry(frame.findObject("bkg_peaking"), f"Peaking background (n_{{peak}} = {n_peak_val:.2f} +/- {n_peak_err:.2f})", "L")
    leg.AddEntry(frame.findObject("bkg_nonpeak"), f"Non-peaking background (n_{{nonpeak}} = {n_nonpeak_val:.2f} +/- {n_nonpeak_err:.2f})", "L")
    r_dummy = ROOT.TLine()
    r_dummy.SetLineWidth(0)
    r_dummy.SetLineColor(0)
    leg.AddEntry(r_dummy, f"n_{{bkg}} = {n_bkg.getVal():.2f} +/- {n_bkg.getError():.2f}", "L")
    leg.AddEntry(r_dummy, f"Ratio r (peaking / total) = {ratio_val:.3f} (fixed)", "L")

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


def plot_toy_fit(toys, model, result, met, mll, mass_point, r=0.088,
                 plotname_prefix="toy_fit", eos_dir="."):
    """
    Plot 1D projections of the toy dataset with the best-fit model overlaid,
    in the same style as fit2D.py.
    """
    n_sig = result.floatParsFinal().find("n_sig")
    n_bkg = result.floatParsFinal().find("n_bkg")

    n_sig_val     = n_sig.getVal()
    n_bkg_val     = n_bkg.getVal()
    n_bkg_err     = n_bkg.getError()
    n_tot_val     = n_sig_val + n_bkg_val
    n_peak_val    = n_bkg_val * r
    n_nonpeak_val = n_bkg_val * (1 - r)
    n_peak_err    = n_bkg_err * r
    n_nonpeak_err = n_bkg_err * (1 - r)

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

        model.plotOn(frame,
                     ROOT.RooFit.Components("sigtot_mll_met_2dpdf"),
                     ROOT.RooFit.Name("signal"),
                     ROOT.RooFit.Normalization(n_tot_val, ROOT.RooAbsReal.NumEvent),
                     ROOT.RooFit.LineColor(ROOT.TColor.GetColor("#bd1f01")),
                     ROOT.RooFit.LineStyle(ROOT.kDashed),
                     ROOT.RooFit.LineWidth(2))

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
                           ratio_val=r, mass_point=mass_point,
                           eos_dir="/eos/user/s/skkwan/www/higgsino/studies/mll-MET-fit-2D/toys")

        make_toy_plot(frame, doLog=False, **plot_kwargs)
        make_toy_plot(frame, doLog=True,  **plot_kwargs)


if __name__ == "__main__" :
    m1 = 650
    m2 = 1
    r = 0.088 # TODO: make this floating

    # Observables
    met = ROOT.RooRealVar("met", "met", 200, 1200)
    mll = ROOT.RooRealVar("m_ll", "m_ll", 60, 120)

    # Generate background toys: get the components as well to keep them in scope in Python
    bkg_expect = 20
    bkg_strength = 1
    num_bkg_toys = bkg_expect * bkg_strength
    bkg_model, bkg_components = get_background_model(r)
    bkg_toys = gen_toys(met, mll, num_bkg_toys, bkg_model)
    print(bkg_toys, type(bkg_toys))
    bkgtoys_tot = bkg_toys.sumEntries()
    print(bkgtoys_tot)

    # Generate signal toys
    sig_expect = 2
    sig_strength = 1
    num_sig_toys = sig_expect * sig_strength
    sig_model, sig_components = get_signal_model(m1=m1, m2=m2)
    sig_toys = gen_toys(met, mll, num_sig_toys, sig_model)

    sigtoys_tot = sig_toys.sumEntries()
    print(sigtoys_tot)

    # Make signal plus background toys, model, and fit
    combined_toys = combine_signal_plus_background_toys(bkg_toys, sig_toys)
    combined_model, combined_components = get_signal_plus_background_model(num_sig_toys, num_bkg_toys, m1, m2, r)
    print(combined_toys)
    comb_result = do_combined_fit(combined_model, combined_toys)

    # Analyze results
    nb_comb_fit = ufloat(comb_result.floatParsFinal().find("n_bkg").getVal(), comb_result.floatParsFinal().find("n_bkg").getError())
    ns_comb_fit = ufloat(comb_result.floatParsFinal().find("n_sig").getVal(), comb_result.floatParsFinal().find("n_sig").getError())

    print(f"{nb_comb_fit:.2f}")
    print(f"{ns_comb_fit:.2f}")

    plot_toy_fit(combined_toys, combined_model, comb_result, met, mll,
                 mass_point=(m1, m2), r=r)
