import ROOT
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

if __name__ == "__main__" :
    m1 = 650
    m2 = 1
    r=0.088

    # Observables
    met = ROOT.RooRealVar("met", "met", 200, 1200)
    mll = ROOT.RooRealVar("m_ll", "m_ll", 60, 120)

    # Generate background toys: get the components as well to keep them in scope in Python
    bkg_expect = 1000
    bkg_strength = 1
    num_bkg_toys = bkg_expect * bkg_strength
    bkg_model, bkg_components = get_background_model(r)
    bkg_toys = gen_toys(met, mll, num_bkg_toys, bkg_model)
    print(bkg_toys, type(bkg_toys))
    bkgtoys_tot = bkg_toys.sumEntries()
    print(bkgtoys_tot)

    # Generate signal toys
    sig_expect = 10
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
