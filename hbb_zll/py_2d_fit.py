# Run in ROOT 6.38 (do not do cmsenv)

import ROOT
from ROOT import RooFit as RF

##### DEFINE FIT OBSERVABLES ####
mll = ROOT.RooRealVar("m_ll", "m_ll", 60, 120) 
met = ROOT.RooRealVar("met", "met", 0, 1200)

# # ##### SIGNAL FIT to signal MC file ######

# #Signal 1d met model
# a_met = ROOT.RooRealVar('a_met', 'a_met', 10, 0, 3000)
# b_met = ROOT.RooRealVar('b_met', 'b_met', 10, 0, 1000)  
# c_met = ROOT.RooRealVar('c_met', 'c_met', 0.5, 0, 1000)  
# e_met = ROOT.RooRealVar('e_met', 'e_met', 1, 0.2, 100) 
# sig_smoid_met = ROOT.RooGenericPdf('sig_smoid_met', '(1-exp(-c_met*met))/(1 + exp((met^e_met-a_met)/b_met))', ROOT.RooArgList(met, a_met, b_met, c_met, e_met))

# #Signal 1d mll model
# mean_mll = ROOT.RooRealVar("mean_mll", "mean_mll", 90, 85, 95)
# sigmal_mll = ROOT.RooRealVar("sigmal_mll", "sigmal_mll", 2, 0.01, 10)
# sigmar_mll = ROOT.RooRealVar("sigmar_mll", "sigmar_mll", 2, 0.01, 10)
# alphal_mll = ROOT.RooRealVar("alphal_mll","alphal_mll", 4, 0.01, 10)
# nl_mll = ROOT.RooRealVar("nl_mll", "nl_mll", 2, 0.01, 100)
# alphar_mll = ROOT.RooRealVar("alphar_mll","alphar_mll", 5, 0.01, 10)
# nr_mll = ROOT.RooRealVar("nr_mll", "nr_mll", 0.01, 0.01, 100)
# sig_dcb_mll = ROOT.RooCrystalBall("sig_dcb_mll", "sig_dcb_mll", mll, mean_mll, sigmal_mll, sigmar_mll, alphal_mll, nl_mll, alphar_mll, nr_mll)

# #Signal 2D model: sigtot_mll_met_2dpdf = sig_smoid_met * sig_dcb_mll
# sigtot_mll_met_2dpdf = ROOT.RooProdPdf("sigtot_dcb_mll_moid_met", "sigtot_dcb_mll_moid_met", [sig_dcb_mll, sig_smoid_met])

# ###### Retrieve signal datasset from signal root file 
# sigfilepath = 'snapshot_TChiZH_650_1_cat_0_batch_0_channel_mm_SR_mll_MET_fit_scheme.root'
# sigfile = ROOT.TFile.Open(sigfilepath, "READ")
# sigtree = sigfile.Get("event_tree")
# weightXyear = ROOT.RooRealVar("weight_nominal_mm", "weight_nominal_mm", -1, 1)
# variables = ROOT.RooArgSet(mll, met, weightXyear)
# sigdataset = ROOT.RooDataSet("sigdataset", "sigdataset", variables, ROOT.RooFit.Import(sigtree), ROOT.RooFit.WeightVar(weightXyear))

# ###### 2D signal fit 
# result = sigtot_mll_met_2dpdf.fitTo(sigdataset, RF.Save()) #where dataset is RooDataSet
# params = result.floatParsFinal()

# print(params)

#### BACKGROUND FIT to data in control region ######

# Background fake mll model in MET dimension
mu_fakemll_met = ROOT.RooRealVar('mu_fakemll_met', 'mu_fakemll_met', 22.2, 10, 30) 
b_fakemll_met = ROOT.RooRealVar('b_fakemll_met', 'b_fakemll_met', 12.2, 5, 30) 
bkgfakemll_met = ROOT.RooGenericPdf("bkgfakemll_met", "bkgfakemll_met", "1/b_fakemll_met * exp(-(@0 - mu_fakemll_met)/b_fakemll_met - exp(-(@0 - mu_fakemll_met)/b_fakemll_met))",
                        ROOT.RooArgList(met, mu_fakemll_met, b_fakemll_met))  
# Background fake mll model in mll dimension: falling exponential (using this for now)
a_fakemll_mll = ROOT.RooRealVar("a_fakemll_mll", "a_fakemll_mll", -0.02, -1, 1) 
bkgfakemll_mll = ROOT.RooExponential("bkgfakemll_mll", "bkgfakemll_mll", mll, a_fakemll_mll)
#Background 2d fakemll model: bkgfakemll_mll_met_2dpdf = bkgfakemll_met * bkgfakemll_mll
bkgfakemll_mll_met_2dpdf = ROOT.RooProdPdf("bkgfakemll_mll_met_2dpdf", "bkgfakemll_mll_met_2dpdf", [bkgfakemll_mll, bkgfakemll_met])


#Background real mll model in met dimension
mu_realmll_met = ROOT.RooRealVar('mu_realmll_met', 'mu_realmll_met', 50, 30, 100)  
b_realmll_met = ROOT.RooRealVar('b_realmll_met', 'b_realmll_met', 29.6, 20, 40)   
bkgrealmll_met = ROOT.RooGenericPdf("bkgrealmll_met", "bkgrealmll_met", "1/b_realmll_met * exp(-(@0 - mu_realmll_met)/b_realmll_met - exp(-(@0 - mu_realmll_met)/b_realmll_met))",
                        ROOT.RooArgList(met, mu_realmll_met, b_realmll_met))  
# Background real mll model in mll dimension: use the same shape as the signal
bkg_mean_mll = ROOT.RooRealVar("bkg_mean_mll", "bkg_mean_mll", 90, 85, 95)
bkg_sigmal_mll = ROOT.RooRealVar("bkg_sigmal_mll", "bkg_sigmal_mll", 2, 0.01, 10)
bkg_sigmar_mll = ROOT.RooRealVar("bkg_sigmar_mll", "bkg_sigmar_mll", 2, 0.01, 10)
bkg_alphal_mll = ROOT.RooRealVar("bkg_alphal_mll", "bkg_alphal_mll", 4, 0.01, 10)
bkg_nl_mll = ROOT.RooRealVar("bkg_nl_mll", "bkg_nl_mll", 2, 0.01, 100)
bkg_alphar_mll = ROOT.RooRealVar("bkg_alphar_mll", "bkg_alphar_mll", 5, 0.01, 10)
bkg_nr_mll = ROOT.RooRealVar("bkg_nr_mll", "bkg_nr_mll", 0.01, 0.01, 100)
bkgrealmll_dcb_mll = ROOT.RooCrystalBall("bkg_dcb_mll", "bkg_dcb_mll", mll, bkg_mean_mll, bkg_sigmal_mll, bkg_sigmar_mll, bkg_alphal_mll, bkg_nl_mll, bkg_alphar_mll, bkg_nr_mll)
#Background 2d realmll model: bkgrealmll_mll_met_2dpdf = bkgrealmll_met * bkgrealmll_dcb_mll
bkgrealmll_mll_met_2dpdf = ROOT.RooProdPdf("bkgrealmll_mll_met_2dpdf", "bkgrealmll_mll_met_2dpdf", [bkgrealmll_dcb_mll, bkgrealmll_met])


#Overall 2D bkg model: bkgtot_mll_met_2dpdf = bkgfakemet_mll_met_2dpdf + ratio_realmll * bkgrealmll_mll_met_2dpdf
ratio_realmll = ROOT.RooRealVar("ratio_realmll", "ratio_realmll", 0.44, 0, 1)
bkgtot_mll_met_2dpdf = ROOT.RooAddPdf("bkgtot_mll_met_2dpdf", "bkgtot_mll_met_2dpdf", [bkgrealmll_mll_met_2dpdf, bkgfakemll_mll_met_2dpdf], [ratio_realmll])

###### Retrive cr data root file ########
crfilepath = 'backgrounds.root'
crfile = ROOT.TFile.Open(crfilepath, "READ")
crtree = crfile.Get("event_tree")
variables = ROOT.RooArgSet(mll, met)
weight = ROOT.RooRealVar("weight_nominal_mm", "weight_nominal_mm", -1, 1)
crdataset = ROOT.RooDataSet("crdataset", "crdataset", variables, ROOT.RooFit.Import(crtree), ROOT.RooFit.WeightVar(weight))

### B only 2D fit to cr root file
result = bkgtot_mll_met_2dpdf.fitTo(crdataset, RF.Save()) #where dataset is RooDataSet
params = result.floatParsFinal()

