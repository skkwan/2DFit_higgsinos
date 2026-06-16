def gen_toys(num_toys, model):

    toys = model.generate(ROOT.RooArgSet(mgg, met), num_toys)

    return toys


def comb_spbtoys(bkgtoys, sigtoys): 
    """
    Combining sig+bkg toys
    """
    spbtoys = bkgtoys.Clone()
    spbtoys.append(sigtoys) 
    return spbtoys

def get_floatsbp_model(nlsp, lsp, sig_index):
    """
    Get the signal plus background model for a specific signal mass point.
    Defines sigtot_dcb_mgg_moid_met_floatspb and bkgtot_mgg_met_floatspb_2dpdf.
    """
    ######### Get signal fitted parameters ##########
    df_fitresult = pd.read_pickle(f"/Users/meraj/bbggAnalysis/mc_ntuples/fitting__ntupleswith__trig_keepevent_photonptthresh/sig_mc/sig_unbinfit_params.pkl")

    #Loading the signal fit parameters, creating model
    df = df_fitresult.loc[sig_index]['df_fitpar']

    if df_fitresult.loc[sig_index]['dataname'] == f'{nlsp}_{lsp}_SR':
        
        #Signal met model
        a_met_value = df.loc[df['Parameter'] == 'a', 'Value'].values[0]
        b_met_value = df.loc[df['Parameter'] == 'b', 'Value'].values[0]
        c_met_value = df.loc[df['Parameter'] == 'c', 'Value'].values[0]
        e_met_value = df.loc[df['Parameter'] == 'e', 'Value'].values[0]
        
        #signal 
        a_floatspb_met = ROOT.RooRealVar(f'a_floatspb_{observablename[1]}', f'a_floatspb_{observablename[1]}', a_met_value)  # Parameter a, with initial value 2 and range [0, 5]
        b_floatspb_met = ROOT.RooRealVar(f'b_floatspb_{observablename[1]}', f'b_floatspb_{observablename[1]}', b_met_value)  # Parameter b, with initial value 1 and range [0.1, 10]
        c_floatspb_met = ROOT.RooRealVar(f'c_floatspb_{observablename[1]}', f'c_floatspb_{observablename[1]}', c_met_value)  # Parameter c, with initial value 0.5 and range [0, 2]
        e_floatspb_met = ROOT.RooRealVar(f'e_floatspb_{observablename[1]}', f'e_floatspb_{observablename[1]}', e_met_value) 
        
        # Create the RooGenericPdf with the specified function
        sig_moid_floatspb_met = ROOT.RooGenericPdf(f'sig_moid_floatspb_{observablename[1]}', '(1-exp(-c_floatspb_met*met))/(1 + exp((met^e_floatspb_met-a_floatspb_met)/b_floatspb_met))', #(1-exp(-c*met^d)) atan(c*met)
                                ROOT.RooArgList(observables[1], a_floatspb_met, b_floatspb_met, c_floatspb_met, e_floatspb_met))
        
        #Signal mgg model
        alphal_mgg_value = df.loc[df['Parameter'] == 'alphal_mgg', 'Value'].values[0]
        alphar_mgg_value = df.loc[df['Parameter'] == 'alphar_mgg', 'Value'].values[0]
        mean_mgg_value = df.loc[df['Parameter'] == 'mean_mgg', 'Value'].values[0]
        nl_mgg_value = df.loc[df['Parameter'] == 'nl_mgg', 'Value'].values[0]
        nr_mgg_value = df.loc[df['Parameter'] == 'nr_mgg', 'Value'].values[0]
        sigmal_mgg_value = df.loc[df['Parameter'] == 'sigmal_mgg', 'Value'].values[0]
        sigmar_mgg_value = df.loc[df['Parameter'] == 'sigmar_mgg', 'Value'].values[0]
        
        mean_floatspb_mgg = ROOT.RooRealVar(f"mean_floatspb_{observablename[0]}", f"mean_floatspb_{observablename[0]}", mean_mgg_value)
        sigmal_floatspb_mgg = ROOT.RooRealVar(f"sigmal_floatspb_{observablename[0]}", f"sigmal_floatspb_{observablename[0]}", sigmal_mgg_value)
        sigmar_floatspb_mgg = ROOT.RooRealVar(f"sigmar_floatspb_{observablename[0]}", f"sigmar_floatspb_{observablename[0]}", sigmar_mgg_value)
        alphal_floatspb_mgg = ROOT.RooRealVar(f"alphal_floatspb_{observablename[0]}",f"alphal_floatspb_{observablename[0]}", alphal_mgg_value)
        nl_floatspb_mgg = ROOT.RooRealVar(f"nl_floatspb_{observablename[0]}", f"nl_floatspb_{observablename[0]}", nl_mgg_value)
        alphar_floatspb_mgg = ROOT.RooRealVar(f"alphar_floatspb_{observablename[0]}",f"alphar_floatspb_{observablename[0]}", alphar_mgg_value)
        nr_floatspb_mgg = ROOT.RooRealVar(f"nr_floatspb_{observablename[0]}", f"nr_floatspb_{observablename[0]}", nr_mgg_value)
        
        sig_dcb_floatspb_mgg = ROOT.RooCrystalBall(f"sig_dcb_floatspb_{observablename[0]}", f"sig_dcb_floatspb_{observablename[0]}", observables[0], mean_floatspb_mgg, sigmal_floatspb_mgg, sigmar_floatspb_mgg, alphal_floatspb_mgg, nl_floatspb_mgg, alphar_floatspb_mgg, nr_floatspb_mgg)
        
        sigtot_mgg_met_floatspb_2dpdf = ROOT.RooProdPdf("sigtot_dcb_mgg_moid_met_floatspb", "sigtot_dcb_mgg_moid_met_floatspb", [sig_dcb_floatspb_mgg, sig_moid_floatspb_met])

        
        #Ptmiss model
        mu_floatspb_fakemet = ROOT.RooRealVar('mu_floatspb_fakemet', 'mu_floatspb_fakemet', 24.4, 10, 30) #24.4, 23, 26
        b_floatspb_fakemet = ROOT.RooRealVar('b_floatspb_fakemet', 'b_floatspb_fakemet', 16.8, 10, 20) #16.8, 15, 17
        bkg_fakemet_gumbel_floatspb_met = ROOT.RooGenericPdf("bkg_fakemet_gumbel_floatspb_met", "bkg_fakemet_gumbel_floatspb_met", "1/b_floatspb_fakemet * exp(-(@0 - mu_floatspb_fakemet)/b_floatspb_fakemet - exp(-(@0 - mu_floatspb_fakemet)/b_floatspb_fakemet))",
                            ROOT.RooArgList(met, mu_floatspb_fakemet, b_floatspb_fakemet))  
        #mgg
        a_fakemet_floatspb_mgg = ROOT.RooRealVar(f"a_fakemet_floatspb_mgg", f"a_fakemet_floatspb_mgg", -0.01,-1,1)
        bkg_exp_fakemet_floatspb_mgg = ROOT.RooExponential(f"bkg_exp_fakemet_floatspb_mgg", f"bkg_exp_fakemet_floatspb_mgg", mgg, a_fakemet_floatspb_mgg)

        bkgfakemet_mgg_met_floatspb_2dpdf = ROOT.RooProdPdf("bkgfakemet_mgg_met_floatspb_2dpdf", "bkgfakemet_mgg_met_floatspb_2dpdf", [bkg_exp_fakemet_floatspb_mgg, bkg_fakemet_gumbel_floatspb_met])
        
        
        #Ptmiss model
        mu_floatspb_realmet = ROOT.RooRealVar('mu_floatspb_realmet', 'mu_floatspb_realmet', 53.4, 30, 70)  #53.4, 45, 62
        b_floatspb_realmet = ROOT.RooRealVar('b_floatspb_realmet', 'b_floatspb_realmet', 29.9, 20, 40)  #29.9, 15, 45
        bkg_realmet_gumbel_floatspb_met = ROOT.RooGenericPdf("bkg_realmet_gumbel_floatspb_met", "bkg_realmet_gumbel_floatspb_met", "1/b_floatspb_realmet * exp(-(@0 - mu_floatspb_realmet)/b_floatspb_realmet - exp(-(@0 - mu_floatspb_realmet)/b_floatspb_realmet))",
                            ROOT.RooArgList(met, mu_floatspb_realmet, b_floatspb_realmet))  
        
        #mgg
        a_realmet_floatspb_mgg = ROOT.RooRealVar(f"a_realmet_floatspb_mgg", f"a_realmet_floatspb_mgg", -0.01,-1,1)
        bkg_exp_realmet_floatspb_mgg = ROOT.RooExponential(f"bkg_exp_realmet_floatspb_mgg", f"bkg_exp_realmet_floatspb_mgg", mgg, a_realmet_floatspb_mgg)

        bkgrealmet_mgg_met_floatspb_2dpdf = ROOT.RooProdPdf("bkgrealmet_mgg_met_floatspb_2dpdf", "bkgrealmet_mgg_met_floatspb_2dpdf", [bkg_exp_realmet_floatspb_mgg, bkg_realmet_gumbel_floatspb_met])

        #Make the overall bkg model
        ratio_floatspb_realmet = ROOT.RooRealVar("ratio_floatspb_realmet", "ratio_floatspb_realmet", 0.1, 0.01, 1)
        bkgtot_mgg_met_floatspb_2dpdf = ROOT.RooAddPdf("bkgtot_mgg_met_floatspb_2dpdf", "bkgtot_mgg_met_floatspb_2dpdf", [bkgrealmet_mgg_met_floatspb_2dpdf, bkgfakemet_mgg_met_floatspb_2dpdf], [ratio_floatspb_realmet])
        n_sig = ROOT.RooRealVar("n_sig", "n_sig", 50, -2, 500)
        n_bkg = ROOT.RooRealVar("n_bkg", "n_bkg", 60, 0, 800)

        sigbkg_floatspb_2dpdf = ROOT.RooAddPdf("sigbkg_floatspb_2dpdf", "sigbkg_floatspb_2dpdf", [sigtot_mgg_met_floatspb_2dpdf, bkgtot_mgg_met_floatspb_2dpdf], [n_sig, n_bkg])
        w_sbp_float = ROOT.RooWorkspace(f"float_2d_sbp{nlsp}_{lsp}")
        w_sbp_float.Import(sigbkg_floatspb_2dpdf)
        return w_sbp_float



def fit_spbtoys(fitmodel, fittoys, signalmass, signalstren):
    """
    Fits the fittoys to fitmodel and analyzes the results.
    """
    result = fitmodel.fitTo(fittoys, RF.Save(), RF.SumW2Error(False))
    cov_np = get_covmatrix(result)
    params = result.floatParsFinal()   #Get plots and gof params
    for j, observable in enumerate(observables):        
    if j == 0:
        ks_statistic_mgg, ks_p_value_mgg = getks(fitmodel, fittoys, observable)
        chi2_ndf_fromplot_mgg, canvaswithfit = plot(fitmodel, fittoys, observable)
        canvaswithfit.Print(f"datafits/{fitmodel.GetName()}_{signalmass}_{signalstren}.pdf[")
    else:
        ks_statistic_met, ks_p_value_met = getks(fitmodel, fittoys, observable)
        chi2_ndf_fromplot_met, canvaswithfit = plot(fitmodel, fittoys, observable)  
    canvaswithfit.Print(f"datafits/{fitmodel.GetName()}_{signalmass}_{signalstren}.pdf")
    
    #Load dataframe with fit results
    df_fitpar = params_to_dataframe(params, ks_statistic_mgg, ks_p_value_mgg, chi2_ndf_fromplot_mgg, ks_statistic_met, ks_p_value_met, chi2_ndf_fromplot_met)
    canvaswithfit.Print(f"datafits/{fitmodel.GetName()}{signalmass}{signalstren}.pdf]")
    return df_fitpar, result



#Study 1
#[nlsp, lsp, signal_index, num_sig (after 20 GeV met cut)]
signal_mass_points = {200: [200,0,0,60],
                      300: [300,0,1,18],
                      500: [500,0,6,4]}
#Analysis table
columns_table = [
    'signal_mass',
    'sig_strength',
   'bkg_strength', 
   'nb_true_tot',
   'nb_fit_tot',
   'nb_true_sr',
   'nb_fit_sr',
   'rel_diff_nb',
   'ns_true_tot', 
    'ns_fit_tot',
   'ns_true_sr',
   'ns_fit_sr',
    'rel_diff_ns'
]
df_table = pd.DataFrame(columns=columns_table) 
# Define 2D window
mgg.setRange('SR', 115, 135)
met.setRange('SR', 50, 500)
#Generate background toys
bkg_expect = 380
bkg_strength = 1
w_sb_fix = get_fix_sb_shape()
sb_shape = w_sb_fix.pdf("bkgtot_mgg_met_sbs_2dpdf")
num_bkg_toys = bkg_strength * bkg_expect

bkg_toys = gen_toys(num_bkg_toys, sb_shape)
bkgtoys_tot = bkg_toys.sumEntries()
bkgtoys_sr = bkg_toys.sumEntries("", 'SR')
nb_true_sr = ufloat(bkgtoys_sr, np.sqrt(bkgtoys_sr))
nb_true_tot = ufloat(bkgtoys_tot, np.sqrt(bkgtoys_tot))
sigstren_list = np.linspace(0, 1, 11)

for signal in signal_mass_points:

    signal_list = signal_mass_points.get(signal)
    sig_expect = signal_list[3]

    for sigstren in sigstren_list:
    #signal toys
    nlsp, lsp, sig_index = signal_list[0], signal_list[1], signal_list[2]
    w_sig = get_signal_shape(nlsp, lsp, sig_index)
    signal_shape = w_sig.pdf("sigtot_dcb_mgg_moid_met")
    num_sig_toys = sig_expect * sigstren
    sig_toys = gen_toys(num_sig_toys, signal_shape)

    sigtoys_tot = sig_toys.sumEntries()
    sigtoys_sr = sig_toys.sumEntries("", 'SR')
    ns_true_sr = ufloat(sigtoys_sr, np.sqrt(sigtoys_sr))
    ns_true_tot = ufloat(sigtoys_tot, np.sqrt(sigtoys_tot))

    #spb toys
    spbtoys = comb_spbtoys(bkg_toys, sig_toys)

    #spb float shape
    w_sbp_float = get_floatsbp_model(nlsp, lsp, sig_index)
    sbp_float_shape = w_sbp_float.pdf("sigbkg_floatspb_2dpdf")

    #fit sbp toys
    df_fitpar, result = fit_spbtoys(sbp_float_shape, spbtoys, signal, int(sigstren*100))
     
    # Compute 2D integrals over the window
    component_pdfs = sbp_float_shape.pdfList()
    for pdf in component_pdfs:
        if pdf.GetName() == 'sigtot_dcb_mgg_moid_met_floatspb':
            sig_int = pdf.createIntegral(ROOT.RooArgSet(mgg, met), ROOT.RooFit.NormSet(ROOT.RooArgSet(mgg, met)), ROOT.RooFit.Range('SR'))
        elif pdf.GetName() == 'bkgtot_mgg_met_floatspb_2dpdf':
            bkg_int = pdf.createIntegral(ROOT.RooArgSet(mgg, met),ROOT.RooFit.NormSet(ROOT.RooArgSet(mgg, met)), ROOT.RooFit.Range('SR'))

    # Get fitted nb and ns events in the window
    n_bkg_window = ROOT.RooProduct("background_yield", "background_yield", ROOT.RooArgList(bkg_int, result.floatParsFinal().find('n_bkg')))
    n_sig_window = ROOT.RooProduct("background_yield", "background_yield", ROOT.RooArgList(sig_int, result.floatParsFinal().find('n_sig')))        
    n_sig_window_err = n_sig_window.getPropagatedError(result)
    n_bkg_window_err = n_bkg_window.getPropagatedError(result)
    ns_fit_sr = ufloat(n_sig_window.getVal(), n_sig_window_err)
    nb_fit_sr = ufloat(n_bkg_window.getVal(), n_bkg_window_err)
    
    # Get nb and ns tot from fit
    nb_fit_tot = ufloat(result.floatParsFinal().find('n_bkg').getVal(), result.floatParsFinal().find('n_bkg').getError())
    ns_fit_tot = ufloat(result.floatParsFinal().find('n_sig').getVal(), result.floatParsFinal().find('n_sig').getError())
    
    rel_diff_ns = (ns_fit_sr - ns_true_sr) / ns_true_sr if ns_true_sr != 0 else -1
    rel_diff_nb = (nb_fit_sr - nb_true_sr) / nb_true_sr if nb_true_sr != 0 else -1
    
    df_table.loc[len(df_table)] = [signal,
                                    sigstren, 
                                   bkg_strength, 
                                   nb_true_tot,
                                   nb_fit_tot,
                                   nb_true_sr,
                                   nb_fit_sr,
                                   rel_diff_nb,
                                   ns_true_tot, 
                                    ns_fit_tot,
                                   ns_true_sr,
                                   ns_fit_sr,
                                    rel_diff_ns]