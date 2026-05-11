import ROOT
from ROOT import RooFit as RF
import numpy as np
from array import array

# Define fit observables
mll = ROOT.RooRealVar("m_ll", "m_ll", 60, 120) 

# Background real mll model in mll dimension: use a simple Gaussian
bkg_mean_mll = ROOT.RooRealVar("bkg_Zpeak_mean_mll", "bkg_Zpeak_mean_mll", 90, 85, 95)
bkg_sigma_mll = ROOT.RooRealVar("bkg_Zpeak_sigma_mll", "bkg_Zpeak_sigma_mll", 2, 0.01, 10)
bkgrealmll_mll = ROOT.RooGaussian("bkg_gaus_mll", "bkg_gaus_mll", mll, bkg_mean_mll, bkg_sigma_mll)


# Test DCB
peak_mean_mll = ROOT.RooRealVar("peak_mean_mll", "peak_mean_mll", 90, 80, 100)
peak_sigmal_mll = ROOT.RooRealVar("peak_sigmal_mll", "peak_sigmal_mll", 5, 1, 20)
peak_sigmar_mll = ROOT.RooRealVar("peak_sigmar_mll", "peak_sigmar_mll", 5, 1, 20)
peak_alphal_mll = ROOT.RooRealVar("peak_alphal_mll", "peak_alphal_mll", 4, 0.01, 10)
peak_nl_mll = ROOT.RooRealVar("peak_nl_mll", "peak_nl_mll", 3, 1, 10)
peak_alphar_mll = ROOT.RooRealVar("peak_alphar_mll","peak_alphar_mll", 5, 0.01, 10)
peak_nr_mll = ROOT.RooRealVar("peak_nr_mll", "peak_nr_mll", 3, 1, 10)
peak_dcb_mll = ROOT.RooCrystalBall("peak_dcb_mll", "peak_dcb_mll", mll, peak_mean_mll, peak_sigmal_mll, peak_sigmar_mll, peak_alphal_mll, peak_nl_mll, peak_alphar_mll, peak_nr_mll)


###### Retrive cr data root file ########
crfilepath = 'backgrounds_CRZ_Zpeak_2018.root'
crfile = ROOT.TFile.Open(crfilepath, "READ")
crtree = crfile.Get("event_tree")
variables = ROOT.RooArgSet(mll)
weight = ROOT.RooRealVar("weight_nominal", "weight_nominal", -1, 1)
crdataset = ROOT.RooDataSet("crdataset", "crdataset", variables, ROOT.RooFit.Import(crtree), ROOT.RooFit.WeightVar(weight))

### B only 2D fit to cr root file
zPeak_CRZ_fit_result = peak_dcb_mll.fitTo(crdataset, RF.Save(), SumW2Error=True) #where dataset is RooDataSet
params = zPeak_CRZ_fit_result.floatParsFinal()
print(params)

f = ROOT.TFile("initial_zPeak_fit_result.root", "RECREATE")
zPeak_CRZ_fit_result.Write("zPeak_CRZ_fit_result")
f.Close()
