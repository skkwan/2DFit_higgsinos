# Run in ROOT 6.38 (do not do cmsenv)

import ROOT
from ROOT import RooFit as RF

##### DEFINE FIT OBSERVABLES ####
mll = ROOT.RooRealVar("m_ll", "m_ll", 60, 120) 
met = ROOT.RooRealVar("met", "met", 0, 1200)

# ##### SIGNAL FIT to signal MC file ######

#Signal 1d met model
a_met = ROOT.RooRealVar('a_met', 'a_met', 10, 0, 3000)
b_met = ROOT.RooRealVar('b_met', 'b_met', 10, 0, 1000)  
c_met = ROOT.RooRealVar('c_met', 'c_met', 0.5, 0, 1000)  
e_met = ROOT.RooRealVar('e_met', 'e_met', 1, 0.2, 100) 
sig_smoid_met = ROOT.RooGenericPdf('sig_smoid_met', '(1-exp(-c_met*met))/(1 + exp((met^e_met-a_met)/b_met))', ROOT.RooArgList(met, a_met, b_met, c_met, e_met))

#Signal 1d mll model
mean_mll = ROOT.RooRealVar("mean_mll", "mean_mll", 90, 85, 95)
sigmal_mll = ROOT.RooRealVar("sigmal_mll", "sigmal_mll", 2, 0.01, 10)
sigmar_mll = ROOT.RooRealVar("sigmar_mll", "sigmar_mll", 2, 0.01, 10)
alphal_mll = ROOT.RooRealVar("alphal_mll","alphal_mll", 4, 0.01, 10)
nl_mll = ROOT.RooRealVar("nl_mll", "nl_mll", 2, 0.01, 100)
alphar_mll = ROOT.RooRealVar("alphar_mll","alphar_mll", 5, 0.01, 10)
nr_mll = ROOT.RooRealVar("nr_mll", "nr_mll", 0.01, 0.01, 100)
sig_dcb_mll = ROOT.RooCrystalBall("sig_dcb_mll", "sig_dcb_mll", mll, mean_mll, sigmal_mll, sigmar_mll, alphal_mll, nl_mll, alphar_mll, nr_mll)

#Signal 2D model: sigtot_mll_met_2dpdf = sig_smoid_met * sig_dcb_mll
sigtot_mll_met_2dpdf = ROOT.RooProdPdf("sigtot_dcb_mll_moid_met", "sigtot_dcb_mll_moid_met", [sig_dcb_mll, sig_smoid_met])

###### Retrieve signal datasset from signal root file 
sigfilepath = 'snapshot_TChiZH_650_1_cat_0_batch_0_channel_mm_SR_mll_MET_fit_scheme.root'
sigfile = ROOT.TFile.Open(sigfilepath, "READ")
sigtree = sigfile.Get("event_tree")
weightXyear = ROOT.RooRealVar("weight_nominal_mm", "weight_nominal_mm", -1, 1)
variables = ROOT.RooArgSet(mll, met, weightXyear)
sigdataset = ROOT.RooDataSet("sigdataset", "sigdataset", variables, ROOT.RooFit.Import(sigtree), ROOT.RooFit.WeightVar(weightXyear))

###### 2D signal fit 
result = sigtot_mll_met_2dpdf.fitTo(sigdataset, RF.Save()) #where dataset is RooDataSet
params = result.floatParsFinal()

print(params)

# #### BACKGROUND FIT to data in control region ######

# #Background fakemet model in met dimension
# mu_fakemet_met = ROOT.RooRealVar('mu_fakemet_met', 'mu_fakemet_met', 22.2, 10, 30) 
# b_fakemet_met = ROOT.RooRealVar('b_fakemet_met', 'b_fakemet_met', 12.2, 5, 30) 
# bkgfakemet_met = ROOT.RooGenericPdf("bkgfakemet_met", "bkgfakemet_met", "1/b_fakemet_met * exp(-(@0 - mu_fakemet_met)/b_fakemet_met - exp(-(@0 - mu_fakemet_met)/b_fakemet_met))",
#                         ROOT.RooArgList(met, mu_fakemet_met, b_fakemet_met))  
# #Background fakemet model in mll dimension
# a_fakemet_mll = ROOT.RooRealVar("a_fakemet_mll", "a_fakemet_mll", -0.02, -1, 1) 
# bkgfakemet_mll = ROOT.RooExponential("bkgfakemet_mll", "bkgfakemet_mll", mll, a_fakemet_mll)
# #Background 2d fakemet model: bkgfakemet_mll_met_2dpdf = bkgfakemet_met * bkgfakemet_mll
# bkgfakemet_mll_met_2dpdf = ROOT.RooProdPdf("bkgfakemet_mll_met_2dpdf", "bkgfakemet_mll_met_2dpdf", [bkgfakemet_mll, bkgfakemet_met])


# #Background realmet model in met dimension
# mu_realmet_met = ROOT.RooRealVar('mu_realmet_met', 'mu_realmet_met', 50, 30, 100)  
# b_realmet_met = ROOT.RooRealVar('b_realmet_met', 'b_realmet_met', 29.6, 20, 40)   
# bkgrealmet_met = ROOT.RooGenericPdf("bkgrealmet_met", "bkgrealmet_met", "1/b_realmet_met * exp(-(@0 - mu_realmet_met)/b_realmet_met - exp(-(@0 - mu_realmet_met)/b_realmet_met))",
#                         ROOT.RooArgList(met, mu_realmet_met, b_realmet_met))  
# #Background realmet model in mll dimension
# a_realmet_mll = ROOT.RooRealVar("a_realmet_mll", "a_realmet_mll", -0.02, -1, 1) 
# bkgrealmet_mll = ROOT.RooExponential("bkgrealmet_mll", "bkgrealmet_mll", mll, a_realmet_mll)
# #Background 2d realmet model: bkgrealmet_mll_met_2dpdf = bkgrealmet_met * bkgrealmet_mll
# bkgrealmet_mll_met_2dpdf = ROOT.RooProdPdf("bkgrealmet_mll_met_2dpdf", "bkgrealmet_mll_met_2dpdf", [bkgrealmet_mll, bkgrealmet_met])


# #Overall 2D bkg model: bkgtot_mll_met_2dpdf = bkgfakemet_mll_met_2dpdf + ratio_realmet * bkgrealmet_mll_met_2dpdf
# ratio_realmet = ROOT.RooRealVar("ratio_realmet", "ratio_realmet", 0.44, 0, 1)
# bkgtot_mll_met_2dpdf = ROOT.RooAddPdf("bkgtot_mll_met_2dpdf", "bkgtot_mll_met_2dpdf", [bkgrealmet_mll_met_2dpdf, bkgfakemet_mll_met_2dpdf], [ratio_realmet])

# ###### Retrive cr data root file ########
# crfilepath = 'CR5_data.root'
# crfile = ROOT.TFile.Open(crfilepath, "READ")
# crtree = crfile.Get("tree")
# variables = ROOT.RooArgSet(mll, met)
# crdataset = ROOT.RooDataSet("crdataset", "crdataset", crtree, variables)

# ### B only 2D fit to cr root file
# result = bkgtot_mll_met_2dpdf.fitTo(crdataset, RF.Save()) #where dataset is RooDataSet
# params = result.floatParsFinal()

