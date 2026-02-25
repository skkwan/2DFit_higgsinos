import ROOT
from ROOT import RooFit as RF

##### DEFINE FIT OBSERVABLES ####
mgg = ROOT.RooRealVar("mgg", "mgg", 100, 200) 
met = ROOT.RooRealVar("met", "met", 0, 400)

##### SIGNAL FIT to signal MC file ######

#Signal 1d met model
a_met = ROOT.RooRealVar('a_met', 'a_met', 10, 0, 3000)
b_met = ROOT.RooRealVar('b_met', 'b_met', 10, 0, 1000)  
c_met = ROOT.RooRealVar('c_met', 'c_met', 0.5, 0, 1000)  
e_met = ROOT.RooRealVar('e_met', 'e_met', 1, 0.2, 100) 
sig_smoid_met = ROOT.RooGenericPdf('sig_smoid_met', '(1-exp(-c_met*met))/(1 + exp((met^e_met-a_met)/b_met))', ROOT.RooArgList(met, a_met, b_met, c_met, e_met))

#Signal 1d mgg model
mean_mgg = ROOT.RooRealVar("mean_mgg", "mean_mgg", 125, 120, 130)
sigmal_mgg = ROOT.RooRealVar("sigmal_mgg", "sigmal_mgg", 2, 0.01, 10)
sigmar_mgg = ROOT.RooRealVar("sigmar_mgg", "sigmar_mgg", 2, 0.01, 10)
alphal_mgg = ROOT.RooRealVar("alphal_mgg","alphal_mgg", 4, 0.01, 10)
nl_mgg = ROOT.RooRealVar("nl_mgg", "nl_mgg", 2, 0.01, 100)
alphar_mgg = ROOT.RooRealVar("alphar_mgg","alphar_mgg", 5, 0.01, 10)
nr_mgg = ROOT.RooRealVar("nr_mgg", "nr_mgg", 0.01, 0.01, 100)
sig_dcb_mgg = ROOT.RooCrystalBall("sig_dcb_mgg", "sig_dcb_mgg", mgg, mean_mgg, sigmal_mgg, sigmar_mgg, alphal_mgg, nl_mgg, alphar_mgg, nr_mgg)

#Signal 2D model: sigtot_mgg_met_2dpdf = sig_smoid_met * sig_dcb_mgg
sigtot_mgg_met_2dpdf = ROOT.RooProdPdf("sigtot_dcb_mgg_moid_met", "sigtot_dcb_mgg_moid_met", [sig_dcb_mgg, sig_smoid_met])

###### Retrieve signal datasset from signal root file 
sigfilepath = 'SR_bbgg_ntuples_SMS-TChiHH_mChi-300_mLSP-0_HToGG_2D_2016.root'
sigfile = ROOT.TFile.Open(sigfilepath, "READ")
sigtree = sigfile.Get("tree")
weightXyear = ROOT.RooRealVar("weightXyear", "weightXyear", -1, 1)
variables = ROOT.RooArgSet(mgg, met, weightXyear)
sigdataset = ROOT.RooDataSet("sigdataset", "sigdataset", sigtree, variables, "", "weightXyear")

###### 2D signal fit 
result = sigtot_mgg_met_2dpdf.fitTo(sigdataset, RF.Save()) #where dataset is RooDataSet
params = result.floatParsFinal()

print(params)

#### BACKGROUND FIT to data in control region ######

#Background fakemet model in met dimension
mu_fakemet_met = ROOT.RooRealVar('mu_fakemet_met', 'mu_fakemet_met', 22.2, 10, 30) 
b_fakemet_met = ROOT.RooRealVar('b_fakemet_met', 'b_fakemet_met', 12.2, 5, 30) 
bkgfakemet_met = ROOT.RooGenericPdf("bkgfakemet_met", "bkgfakemet_met", "1/b_fakemet_met * exp(-(@0 - mu_fakemet_met)/b_fakemet_met - exp(-(@0 - mu_fakemet_met)/b_fakemet_met))",
                        ROOT.RooArgList(met, mu_fakemet_met, b_fakemet_met))  
#Background fakemet model in mgg dimension
a_fakemet_mgg = ROOT.RooRealVar("a_fakemet_mgg", "a_fakemet_mgg", -0.02, -1, 1) 
bkgfakemet_mgg = ROOT.RooExponential("bkgfakemet_mgg", "bkgfakemet_mgg", mgg, a_fakemet_mgg)
#Background 2d fakemet model: bkgfakemet_mgg_met_2dpdf = bkgfakemet_met * bkgfakemet_mgg
bkgfakemet_mgg_met_2dpdf = ROOT.RooProdPdf("bkgfakemet_mgg_met_2dpdf", "bkgfakemet_mgg_met_2dpdf", [bkgfakemet_mgg, bkgfakemet_met])


#Background realmet model in met dimension
mu_realmet_met = ROOT.RooRealVar('mu_realmet_met', 'mu_realmet_met', 50, 30, 100)  
b_realmet_met = ROOT.RooRealVar('b_realmet_met', 'b_realmet_met', 29.6, 20, 40)   
bkgrealmet_met = ROOT.RooGenericPdf("bkgrealmet_met", "bkgrealmet_met", "1/b_realmet_met * exp(-(@0 - mu_realmet_met)/b_realmet_met - exp(-(@0 - mu_realmet_met)/b_realmet_met))",
                        ROOT.RooArgList(met, mu_realmet_met, b_realmet_met))  
#Background realmet model in mgg dimension
a_realmet_mgg = ROOT.RooRealVar("a_realmet_mgg", "a_realmet_mgg", -0.02, -1, 1) 
bkgrealmet_mgg = ROOT.RooExponential("bkgrealmet_mgg", "bkgrealmet_mgg", mgg, a_realmet_mgg)
#Background 2d realmet model: bkgrealmet_mgg_met_2dpdf = bkgrealmet_met * bkgrealmet_mgg
bkgrealmet_mgg_met_2dpdf = ROOT.RooProdPdf("bkgrealmet_mgg_met_2dpdf", "bkgrealmet_mgg_met_2dpdf", [bkgrealmet_mgg, bkgrealmet_met])


#Overall 2D bkg model: bkgtot_mgg_met_2dpdf = bkgfakemet_mgg_met_2dpdf + ratio_realmet * bkgrealmet_mgg_met_2dpdf
ratio_realmet = ROOT.RooRealVar("ratio_realmet", "ratio_realmet", 0.44, 0, 1)
bkgtot_mgg_met_2dpdf = ROOT.RooAddPdf("bkgtot_mgg_met_2dpdf", "bkgtot_mgg_met_2dpdf", [bkgrealmet_mgg_met_2dpdf, bkgfakemet_mgg_met_2dpdf], [ratio_realmet])

###### Retrive cr data root file ########
crfilepath = 'CR5_data.root'
crfile = ROOT.TFile.Open(crfilepath, "READ")
crtree = crfile.Get("tree")
variables = ROOT.RooArgSet(mgg, met)
crdataset = ROOT.RooDataSet("crdataset", "crdataset", crtree, variables)

### B only 2D fit to cr root file
result = bkgtot_mgg_met_2dpdf.fitTo(crdataset, RF.Save()) #where dataset is RooDataSet
params = result.floatParsFinal()

