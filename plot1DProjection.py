import os, ROOT
import cmsstyle as CMS


def addOverflow(h: ROOT.TH1F) -> ROOT.TH1F:
    """
    Add overflow to a histogram
    """
    h.SetBinContent(h.GetNbinsX(), h.GetBinContent(h.GetNbinsX()) + h.GetBinContent(h.GetNbinsX() + 1))
    return h

def normalizeHist(h: ROOT.TH1F) -> None:
    """
    Normalize a histogram
    """
    h.Scale(1/h.Integral())

def getTChainRDF(listOfFiles: list[str], treeName: str) -> tuple[ROOT.TChain, ROOT.RDataFrame]:
    """
    Return a TChain and RDataFrame from a list of file paths
    """
    ch = ROOT.TChain(treeName)
    for file in listOfFiles:
        ch.Add(file)
    df = ROOT.RDataFrame(ch)
    return ch, df

block1 = [
	"SR_bbgg_ntuples_SMS-TChiHH_mChi-300_mLSP-0_HToGG_2D_2016.root"
]

block2 = [
    "CR5_data.root"
]

ch1, df1 = getTChainRDF(block1, "tree")
ch2, df2 = getTChainRDF(block2, "tree")
# df = df.Define("gen_deltaPhi_ll_ptmiss", "compute_deltaPhi(gen_leps_p4.Phi(), gen_p4_ptmiss.Phi())")

mgg = ROOT.RooRealVar("mgg", "mgg", 60, 260) 
met = ROOT.RooRealVar("met", "met", 0, 400)
weightXyear = ROOT.RooRealVar("weightXyear", "weightXyear", -1, 1)

variablesInfo = [
    ["mgg", "m_{#gamma#gamma} / GeV", 40, 60., 260., mgg],
    ["met", "MET / GeV", 40, 0., 400., met], 
    # ["weightXyear", "WeightXyear", 40, 0., 10.],
]

###### Retrieve signal datasset from signal root file 
sigfilepath = 'SR_bbgg_ntuples_SMS-TChiHH_mChi-300_mLSP-0_HToGG_2D_2016.root'
sigfile = ROOT.TFile.Open(sigfilepath, "READ")
sigtree = sigfile.Get("tree")
variables = ROOT.RooArgSet(mgg, met, weightXyear)
sigdataset = ROOT.RooDataSet("sigdataset", "sigdataset", sigtree, variables, "", "weightXyear")

###### Retrive cr data root file ########
crfilepath = 'CR5_data.root'
crfile = ROOT.TFile.Open(crfilepath, "READ")
crtree = crfile.Get("tree")
variables = ROOT.RooArgSet(mgg, met)
crdataset = ROOT.RooDataSet("crdataset", "crdataset", crtree, variables)


d = {}
dPdf = {}

colors = [ROOT.TColor.GetColor("#832db6"),  # purple
          ROOT.TColor.GetColor("#bd1f01"),  # red
          ROOT.TColor.GetColor("#717581")]  # grey

d["signal"] = {}
d["signal"]["dataframe"] = df1
d["signal"]["color"] = ROOT.TColor.GetColor("#5790fc") # blue 
d["signal"]["label"] = "Signal"
d["signal"]["linewidth"] = 2
d["signal"]["dataset"] = sigdataset
d["signal"]["name"] = "signal"
print("Done getting first chain")

d["bkg"] = {}
d["bkg"]["dataframe"] = df2
d["bkg"]["color"] = ROOT.TColor.GetColor("#9c9ca1") # grey
d["bkg"]["label"] = "Background"
d["bkg"]["linewidth"] = 2
d["bkg"]["dataset"] = crdataset 
d["bkg"]["name"] = "background"

# a_met     = 632.925      +/-  2.80634   (limited)
# alphal_mgg        = 1.66588      +/-  0.00851659        (limited)
# alphar_mgg        = 1.17512      +/-  0.00735136        (limited)
# b_met     = 18.2921      +/-  0.307231  (limited)
# c_met     = 976.268      +/-  625.824   (limited)
# e_met     = 0.893406     +/-  0.0190059 (limited)
# mean_mgg          = 124.943      +/-  0.0114606 (limited)
# nl_mgg    = 2.5573       +/-  0.0361191 (limited)
# nr_mgg    = 32.2549      +/-  0.107155  (limited)
# sigmal_mgg        = 1.52331      +/-  0.00821082        (limited)
# sigmar_mgg        = 1.18388      +/-  0.00737473        (limited)
mean_mgg = ROOT.RooRealVar("mean_mgg", "mean_mgg", 125, 120, 130)
sigmal_mgg = ROOT.RooRealVar("sigmal_mgg", "sigmal_mgg", 2, 0.01, 10)
sigmar_mgg = ROOT.RooRealVar("sigmar_mgg", "sigmar_mgg", 2, 0.01, 10)
alphal_mgg = ROOT.RooRealVar("alphal_mgg","alphal_mgg", 4, 0.01, 10)
nl_mgg = ROOT.RooRealVar("nl_mgg", "nl_mgg", 2, 0.01, 100)
alphar_mgg = ROOT.RooRealVar("alphar_mgg","alphar_mgg", 5, 0.01, 10)
nr_mgg = ROOT.RooRealVar("nr_mgg", "nr_mgg", 0.01, 0.01, 100)

a_met = ROOT.RooRealVar('a_met', 'a_met', 10, 0, 3000)
b_met = ROOT.RooRealVar('b_met', 'b_met', 10, 0, 1000)  
c_met = ROOT.RooRealVar('c_met', 'c_met', 0.5, 0, 1000)  
e_met = ROOT.RooRealVar('e_met', 'e_met', 1, 0.2, 100) 
a_met.setVal(632.925)
b_met.setVal(18.2921)
c_met.setVal(976.268)
e_met.setVal(0.893406)
sig_smoid_met = ROOT.RooGenericPdf('sig_smoid_met', '(1-exp(-c_met*met))/(1 + exp((met^e_met-a_met)/b_met))', ROOT.RooArgList(met, a_met, b_met, c_met, e_met))

mean_mgg.setVal(124.943)
sigmal_mgg.setVal(1.52331)
sigmar_mgg.setVal(1.18388)
alphal_mgg.setVal(1.66588)
alphar_mgg.setVal(1.17512)
nl_mgg.setVal(2.5573)
alphar_mgg.setVal(1.17512)
nr_mgg.setVal(32.2549)
sig_dcb_mgg = ROOT.RooCrystalBall("sig_dcb_mgg", "sig_dcb_mgg", mgg, mean_mgg, sigmal_mgg, sigmar_mgg, alphal_mgg, nl_mgg, alphar_mgg, nr_mgg)

#Signal 2D model: sigtot_mgg_met_2dpdf = sig_smoid_met * sig_dcb_mgg
sigtot_mgg_met_2dpdf = ROOT.RooProdPdf("sigtot_dcb_mgg_moid_met", "sigtot_dcb_mgg_moid_met", [sig_dcb_mgg, sig_smoid_met])

# Fill dPdf by hand 
dPdf["signal"] = {}
for varInfo in variablesInfo:
    dPdf["signal"][f"hist_{varInfo[0]}"] = {}
dPdf["signal"]["hist_mgg"]["var"] = mgg 
dPdf["signal"]["hist_mgg"]["pdf"] = sigtot_mgg_met_2dpdf
dPdf["signal"]["hist_mgg"]["color"] = ROOT.TColor.GetColor("#f89c20") # light orange
dPdf["signal"]["hist_mgg"]["label"] = "Signal 2D fit: 1D projection"
dPdf["signal"]["hist_mgg"]["name"] = "signal1Dproj_mgg"
dPdf["signal"]["hist_mgg"]["linestyle"] = 1
dPdf["signal"]["hist_mgg"]["ymax"] = 10

dPdf["signal"]["hist_met"]["var"] = met 
dPdf["signal"]["hist_met"]["pdf"] = sigtot_mgg_met_2dpdf
dPdf["signal"]["hist_met"]["color"] = ROOT.TColor.GetColor("#f89c20") # light orange
dPdf["signal"]["hist_met"]["label"] = "Signal 2D fit: 1D projection"
dPdf["signal"]["hist_met"]["name"] = "signal1Dproj_met"
dPdf["signal"]["hist_met"]["linestyle"] = 1
dPdf["signal"]["hist_met"]["ymax"] = 0.8

#Background fakemet model in met dimension
mu_fakemet_met = ROOT.RooRealVar('mu_fakemet_met', 'mu_fakemet_met', 22.2, 10, 30) 
mu_fakemet_met.setVal(22.1267)
b_fakemet_met = ROOT.RooRealVar('b_fakemet_met', 'b_fakemet_met', 12.2, 5, 30) 
b_fakemet_met.setVal(12.2339)
bkgfakemet_met = ROOT.RooGenericPdf("bkgfakemet_met", "bkgfakemet_met", "1/b_fakemet_met * exp(-(@0 - mu_fakemet_met)/b_fakemet_met - exp(-(@0 - mu_fakemet_met)/b_fakemet_met))",
                        ROOT.RooArgList(met, mu_fakemet_met, b_fakemet_met))  
#Background fakemet model in mgg dimension
a_fakemet_mgg = ROOT.RooRealVar("a_fakemet_mgg", "a_fakemet_mgg", -0.02, -1, 1) 
a_fakemet_mgg.setVal(-0.020442)
bkgfakemet_mgg = ROOT.RooExponential("bkgfakemet_mgg", "bkgfakemet_mgg", mgg, a_fakemet_mgg)
#Background 2d fakemet model: bkgfakemet_mgg_met_2dpdf = bkgfakemet_met * bkgfakemet_mgg
bkgfakemet_mgg_met_2dpdf = ROOT.RooProdPdf("bkgfakemet_mgg_met_2dpdf", "bkgfakemet_mgg_met_2dpdf", [bkgfakemet_mgg, bkgfakemet_met])

#Background realmet model in met dimension
mu_realmet_met = ROOT.RooRealVar('mu_realmet_met', 'mu_realmet_met', 50, 30, 100)  
mu_realmet_met.setVal(49.4322)
b_realmet_met = ROOT.RooRealVar('b_realmet_met', 'b_realmet_met', 29.6, 20, 40)   
b_realmet_met.setVal(29.5402)

bkgrealmet_met = ROOT.RooGenericPdf("bkgrealmet_met", "bkgrealmet_met", "1/b_realmet_met * exp(-(@0 - mu_realmet_met)/b_realmet_met - exp(-(@0 - mu_realmet_met)/b_realmet_met))",
                        ROOT.RooArgList(met, mu_realmet_met, b_realmet_met))  
#Background realmet model in mgg dimension
a_realmet_mgg = ROOT.RooRealVar("a_realmet_mgg", "a_realmet_mgg", -0.02, -1, 1) 
a_realmet_mgg.setVal(-0.0206688)

bkgrealmet_mgg = ROOT.RooExponential("bkgrealmet_mgg", "bkgrealmet_mgg", mgg, a_realmet_mgg)
#Background 2d realmet model: bkgrealmet_mgg_met_2dpdf = bkgrealmet_met * bkgrealmet_mgg
bkgrealmet_mgg_met_2dpdf = ROOT.RooProdPdf("bkgrealmet_mgg_met_2dpdf", "bkgrealmet_mgg_met_2dpdf", [bkgrealmet_mgg, bkgrealmet_met])

#Overall 2D bkg model: bkgtot_mgg_met_2dpdf = bkgfakemet_mgg_met_2dpdf + ratio_realmet * bkgrealmet_mgg_met_2dpdf
ratio_realmet = ROOT.RooRealVar("ratio_realmet", "ratio_realmet", 0.44, 0, 1)
ratio_realmet.setVal(0.374335)
bkgtot_mgg_met_2dpdf = ROOT.RooAddPdf("bkgtot_mgg_met_2dpdf", "bkgtot_mgg_met_2dpdf", [bkgrealmet_mgg_met_2dpdf, bkgfakemet_mgg_met_2dpdf], [ratio_realmet])

dPdf["bkg"] = {}
for varInfo in variablesInfo:
    dPdf["bkg"][f"hist_{varInfo[0]}"] = {} 
dPdf["bkg"]["hist_mgg"]["var"] = mgg
dPdf["bkg"]["hist_mgg"]["pdf"] = bkgtot_mgg_met_2dpdf
dPdf["bkg"]["hist_mgg"]["color"] = ROOT.TColor.GetColor("#964a8b") # purple
dPdf["bkg"]["hist_mgg"]["label"] = "Background 2D fit: 1D proj"
dPdf["bkg"]["hist_mgg"]["name"] = "bkg1Dproj"
dPdf["bkg"]["hist_mgg"]["linestyle"] = 1
dPdf["bkg"]["hist_mgg"]["ymax"] = 40


dPdf["bkg"]["hist_met"]["var"] = met
dPdf["bkg"]["hist_met"]["pdf"] = bkgtot_mgg_met_2dpdf
dPdf["bkg"]["hist_met"]["color"] = ROOT.TColor.GetColor("#964a8b") # purple
dPdf["bkg"]["hist_met"]["label"] = "Background 2D fit: 1D proj"
dPdf["bkg"]["hist_met"]["name"] = "bkg1Dproj"
dPdf["bkg"]["hist_met"]["linestyle"] = 1
dPdf["bkg"]["hist_met"]["ymax"] = 100

for varInfo in variablesInfo:
    variable = varInfo[0]
    xlabel = varInfo[1]
    nBins = varInfo[2]
    xmin = varInfo[3]
    xmax = varInfo[4]
    rooVar = varInfo[5]

    legYmin = 0.89 - 0.05 * 4
    legYmax = 0.89
    legXmin = 0.5
    legXmax = 0.5 + 0.4
    if ("eta" in variable) or ("phi" in variable):
        legYmin = 0.40 - 0.05 * 4
        legYmax = 0.40
        legXmin = 0.3
        legXmax = 0.4 + 0.4


    # Plot both signal and background points
    for key in d:
        # Each variable only needs one frame
        frame = rooVar.frame(40)

        print("Doing", key)
        leg = CMS.cmsLeg(legXmin, legYmin, legXmax, legYmax, textSize=0.03)
        leg.SetHeader("bb#gamma#gamma: results of 2D fit")

        CMS.SetLumi("")

        # print(y_maxima)
        canv = CMS.cmsCanvas(variable, xmin, xmax, 0, dPdf[key][f"hist_{variable}"]["ymax"], # max(y_maxima) * 1.25,
                            nameXaxis = xlabel,
                            nameYaxis = 'Shape (A.U.)',
                            square = CMS.kSquare, extraSpace=0.05, yTitOffset=1.6, iPos=0)
        canv.SetRightMargin(0.05)
        # CMS.SetExtraText("Private work (CMS simulation)")
        # CMS.SetCmsTextFont(52)
        # CMS.SetCmsTextSize(0.75*0.76)
        CMS.UpdatePad(canv)

        y_maxima = []

        # Plot the dataset fitted

        d[key]["dataset"].plotOn(frame, ROOT.RooFit.Name(d[key]["name"]),
                                    ROOT.RooFit.LineColor(d[key]["color"]),
                                    ROOT.RooFit.LineWidth(d[key]["linewidth"]),
                                    ROOT.RooFit.MarkerColor(d[key]["color"]),
                                    ROOT.RooFit.MarkerSize(1))
        leg.AddEntry(frame.findObject(d[key]["name"]), d[key]["label"])

        # Plot the PDF with the fitted parameters
        # See all draw options: https://root.cern.ch/doc/v636/classRooAbsPdf.html#aa0f2f98d89525302a06a1b7f1b0c2aa6 
        dPdf[key][f"hist_{variable}"]["pdf"].plotOn(frame,
                                                    ROOT.RooFit.Name(dPdf[key][f"hist_{variable}"]["name"]),
                                                    ROOT.RooFit.LineColor(dPdf[key][f"hist_{variable}"]["color"]),
                                                    ROOT.RooFit.LineWidth(2),
                                                    ROOT.RooFit.LineStyle(dPdf[key][f"hist_{variable}"]["linestyle"]),
                                                    ROOT.RooFit.MarkerSize(0))
        leg.AddEntry(frame.findObject(dPdf[key][f"hist_{variable}"]["name"]), dPdf[key][f"hist_{variable}"]["label"])
        frame.Draw("SAME")

        CMS.cmsObjectDraw(leg)

        CMS.UpdatePad(canv)

        outdir = f"/eos/user/s/skkwan/www/higgsino/studies/from-meraj-fit-2D"
        sampleName = d[key]["name"]
        canv.SaveAs(f"{sampleName}-{variable}.pdf")
        canv.SaveAs(f"{sampleName}-{variable}.png")
        os.system(f"mv {sampleName}-{variable}.* {outdir}")

        del(canv)
        del(leg)

