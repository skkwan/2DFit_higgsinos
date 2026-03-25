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
    df = ROOT.RDataFrame(ch, {"m_ll", "met", "weight_mm_nominal"})
    return ch, df

block1 = [
	"snapshot_TChiZH_650_1_cat_0_batch_0_channel_mm_SR_mll_MET_fit_scheme.root"
]

block2 = [
    "backgrounds.root"
]

ch1, df1 = getTChainRDF(block1, "event_tree")
ch2, df2 = getTChainRDF(block2, "event_tree")
# df = df.Define("gen_deltaPhi_ll_ptmiss", "compute_deltaPhi(gen_leps_p4.Phi(), gen_p4_ptmiss.Phi())")

mll = ROOT.RooRealVar("m_ll", "m_ll", 60, 120) 
met = ROOT.RooRealVar("met", "met", 0, 1200)
weightXyear = ROOT.RooRealVar("weight_nominal_mm", "weight_nominal_mm", -1, 1)

variablesInfo = [
    ["mll", "m(ll) / GeV", 40, 60., 120., mll],
    ["met", "MET / GeV", 120, 0., 1200., met], 
    # ["weightXyear", "WeightXyear", 40, 0., 10.],
]

##### Retrieve signal datasset from signal root file 
sigfilepath = 'snapshot_TChiZH_650_1_cat_0_batch_0_channel_mm_SR_mll_MET_fit_scheme.root'
sigfile = ROOT.TFile.Open(sigfilepath, "READ")
sigtree = sigfile.Get("event_tree")
variables = ROOT.RooArgSet(mll, met, weightXyear)
sigdataset = ROOT.RooDataSet("sigdataset", "sigdataset", variables, ROOT.RooFit.Import(sigtree), ROOT.RooFit.WeightVar(weightXyear))

###### Retrive cr data root file ########
crfilepath = 'backgrounds.root'
crfile = ROOT.TFile.Open(crfilepath, "READ")
crtree = crfile.Get("event_tree")
variables = ROOT.RooArgSet(mll, met)
bkgdataset = ROOT.RooDataSet("bkgdataset", "bkgdataset", variables, ROOT.RooFit.Import(crtree), ROOT.RooFit.WeightVar(weightXyear))


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
d["bkg"]["dataset"] = bkgdataset 
d["bkg"]["name"] = "background"

# a_met     = 565.321      +/-  1759.26   (limited)
# alphal_mll        = 1.21829      +/-  6.15304   (limited)
# alphar_mll        = 1.16314      +/-  7.52947   (limited)
# b_met     = 49.9658      +/-  148.592   (limited)
# c_met     = 1.19365e-05  +/-  0.0199096 (limited)
# e_met     = 0.988535     +/-  0.371576  (limited)
# mean_mll          = 90.6656      +/-  7.09258   (limited)
# nl_mll    = 2.35992      +/-  21.0617   (limited)
# nr_mll    = 3.21394      +/-  89.8507   (limited)
# sigmal_mll        = 2.20338      +/-  6.86993   (limited)
# sigmar_mll        = 2.40557      +/-  5.56478   (limited)

## OLD: sigmoid for signal met
a_met = ROOT.RooRealVar('a_met', 'a_met', 10, 0, 3000)
b_met = ROOT.RooRealVar('b_met', 'b_met', 10, 0, 1000)  
c_met = ROOT.RooRealVar('c_met', 'c_met', 0.5, 0, 1000)  
e_met = ROOT.RooRealVar('e_met', 'e_met', 1, 0.2, 100) 
a_met.setVal(565.321)
b_met.setVal(49.9658)
c_met.setVal(1.19365e-05)
e_met.setVal(0.988535)
sig_smoid_met = ROOT.RooGenericPdf('sig_smoid_met', '(1-exp(-c_met*met))/(1 + exp((met^e_met-a_met)/b_met))', ROOT.RooArgList(met, a_met, b_met, c_met, e_met))

### NEW TEST: DCB for signal met
# mean_met = ROOT.RooRealVar("mean_met", "mean_met", 400, 0, 1200)
# sigmal_met = ROOT.RooRealVar("sigmal_met", "sigmal_met", 2, 0.01, 10)
# sigmar_met = ROOT.RooRealVar("sigmar_met", "sigmar_met", 2, 0.01, 10)
# alphal_met = ROOT.RooRealVar("alphal_met","alphal_met", 4, 0.01, 10)
# nl_met = ROOT.RooRealVar("nl_met", "nl_met", 100, 10, 200)
# alphar_met = ROOT.RooRealVar("alphar_met","alphar_met", 5, 0.01, 10)
# nr_met = ROOT.RooRealVar("nr_met", "nr_met", 100, 10, 200)
# mean_met.setVal(379)
# sigmal_met.setVal(6.5001)
# sigmar_met.setVal(4.91442)
# alphal_met.setVal(0.105149)
# alphar_met.setVal(0.0342216)
# nl_met.setVal(11.1361)
# nr_met.setVal(38.5208)
# sig_dcb_met = ROOT.RooCrystalBall("sig_dcb_met", "sig_dcb_met", met, mean_met, sigmal_met, sigmar_met, alphal_met, nl_met, alphar_met, nr_met)


mean_mll = ROOT.RooRealVar("mean_mll", "mean_mll", 90, 85, 95)
sigmal_mll = ROOT.RooRealVar("sigmal_mll", "sigmal_mll", 2, 0.01, 10)
sigmar_mll = ROOT.RooRealVar("sigmar_mll", "sigmar_mll", 2, 0.01, 10)
alphal_mll = ROOT.RooRealVar("alphal_mll","alphal_mll", 4, 0.01, 10)
nl_mll = ROOT.RooRealVar("nl_mll", "nl_mll", 2, 0.01, 100)
alphar_mll = ROOT.RooRealVar("alphar_mll","alphar_mll", 5, 0.01, 10)
nr_mll = ROOT.RooRealVar("nr_mll", "nr_mll", 0.01, 0.01, 100)

mean_mll.setVal(90.6656)
sigmal_mll.setVal(2.20339)
sigmar_mll.setVal(2.40555)
alphal_mll.setVal(1.21829)
alphar_mll.setVal(1.16313)
nl_mll.setVal(2.35992)
nr_mll.setVal(3.21394)
sig_dcb_mll = ROOT.RooCrystalBall("sig_dcb_mll", "sig_dcb_mll", mll, mean_mll, sigmal_mll, sigmar_mll, alphal_mll, nl_mll, alphar_mll, nr_mll)

#Signal 2D model: sigtot_mll_met_2dpdf = sig_smoid_met * sig_dcb_mll
# TODO: testing DCB for signal met instead of sigmoid
# sigtot_mll_met_2dpdf = ROOT.RooProdPdf("sigtot_dcb_mll_moid_met", "sigtot_dcb_mll_moid_met", [sig_dcb_mll, sig_smoid_met])
sigtot_mll_met_2dpdf = ROOT.RooProdPdf("sigtot_dcb_mll_dcb_met", "sigtot_dcb_mll_dcb_met", [sig_dcb_mll, sig_smoid_met])

# Fill dPdf by hand 
dPdf["signal"] = {}
for varInfo in variablesInfo:
    dPdf["signal"][f"hist_{varInfo[0]}"] = {}
dPdf["signal"]["hist_mll"]["var"] = mll 
dPdf["signal"]["hist_mll"]["pdf"] = sigtot_mll_met_2dpdf
dPdf["signal"]["hist_mll"]["color"] = ROOT.TColor.GetColor("#f89c20") # light orange
dPdf["signal"]["hist_mll"]["label"] = "Signal 2D fit: 1D projection"
dPdf["signal"]["hist_mll"]["name"] = "signal1Dproj_mll"
dPdf["signal"]["hist_mll"]["linestyle"] = 1
dPdf["signal"]["hist_mll"]["ymax"] = 0.6

dPdf["signal"]["hist_met"]["var"] = met 
dPdf["signal"]["hist_met"]["pdf"] = sigtot_mll_met_2dpdf
dPdf["signal"]["hist_met"]["color"] = ROOT.TColor.GetColor("#f89c20") # light orange
dPdf["signal"]["hist_met"]["label"] = "Signal 2D fit: 1D projection"
dPdf["signal"]["hist_met"]["name"] = "signal1Dproj_met"
dPdf["signal"]["hist_met"]["linestyle"] = 1
dPdf["signal"]["hist_met"]["ymax"] = 0.4


# a_fakemll_mll     = 0.755871     +/-  1.30419   (limited)
# b_fakemll_met     = 28.9782      +/-  24.4727   (limited)
# b_realmll_met     = 40   +/-  0.480795  (limited)
# bkg_alphal_mll    = 1.67984      +/-  0.997158  (limited)
# bkg_alphar_mll    = 0.3251       +/-  0.687871  (limited)
# bkg_mean_mll      = 91.3098      +/-  0.0396237 (limited)
# bkg_nl_mll        = 0.0694016    +/-  0.123819  (limited)
# bkg_nr_mll        = 0.387898     +/-  0.358384  (limited)
# bkg_sigmal_mll    = 0.0100032    +/-  0.421628  (limited)
# bkg_sigmar_mll    = 0.0108922    +/-  0.174022  (limited)
# mu_fakemll_met    = 29.8716      +/-  19.3964   (limited)
# mu_realmll_met    = 100  +/-  1.59049   (limited)
# ratio_realmll     = 1    +/-  0.0433235 (limited)


#Background real mll model in met dimension
mu_fakemll_met = ROOT.RooRealVar('mu_fakemll_met', 'mu_fakemll_met', 210, 10, 400) 
mu_fakemll_met.setVal(242.523) # by-hand 225
b_fakemll_met = ROOT.RooRealVar('b_fakemll_met', 'b_fakemll_met', 12.2, 5, 60) 
b_fakemll_met.setVal(36.6516) # by-hand 40 
bkgfakemll_met = ROOT.RooGenericPdf("bkgfakemll_met", "bkgfakemll_met", "1/b_fakemll_met * exp(-(@0 - mu_fakemll_met)/b_fakemll_met - exp(-(@0 - mu_fakemll_met)/b_fakemll_met))",
                        ROOT.RooArgList(met, mu_fakemll_met, b_fakemll_met))  
#Background fakemll model in mll dimension
a_fakemll_mll = ROOT.RooRealVar("a_fakemll_mll", "a_fakemll_mll", -0.03, -1, 1) 
a_fakemll_mll.setVal(-0.0210225) # -0.03 by hand
bkgfakemll_mll = ROOT.RooExponential("bkgfakemll_mll", "bkgfakemll_mll", mll, a_fakemll_mll)
#Background 2d fakemll model: bkgfakemll_mll_met_2dpdf = bkgfakemll_met * bkgfakemll_mll
bkgfakemll_mll_met_2dpdf = ROOT.RooProdPdf("bkgfakemll_mll_met_2dpdf", "bkgfakemll_mll_met_2dpdf", [bkgfakemll_mll, bkgfakemll_met])

#Background realmll model in met dimension
mu_realmll_met = ROOT.RooRealVar('mu_realmll_met', 'mu_realmll_met', 50, 30, 400)  
mu_realmll_met.setVal(246.411)
b_realmll_met = ROOT.RooRealVar('b_realmll_met', 'b_realmll_met', 29.6, 20, 100)   
b_realmll_met.setVal(41.1751)
bkgrealmll_met = ROOT.RooGenericPdf("bkgrealmll_met", "bkgrealmll_met", "1/b_realmll_met * exp(-(@0 - mu_realmll_met)/b_realmll_met - exp(-(@0 - mu_realmll_met)/b_realmll_met))",
                        ROOT.RooArgList(met, mu_realmll_met, b_realmll_met))  

# Background realmll in mll dimension
bkg_mean_mll = ROOT.RooRealVar("bkg_mean_mll", "bkg_mean_mll", 90, 85, 95)
bkg_sigmal_mll = ROOT.RooRealVar("bkg_sigmal_mll", "bkg_sigmal_mll", 2, 0.01, 10)
bkg_sigmar_mll = ROOT.RooRealVar("bkg_sigmar_mll", "bkg_sigmar_mll", 2, 0.01, 10)
bkg_alphal_mll = ROOT.RooRealVar("bkg_alphal_mll", "bkg_alphal_mll", 4, 0.01, 10)
bkg_nl_mll = ROOT.RooRealVar("bkg_nl_mll", "bkg_nl_mll", 2, 0.01, 100)
bkg_alphar_mll = ROOT.RooRealVar("bkg_alphar_mll", "bkg_alphar_mll", 5, 0.01, 10)
bkg_nr_mll = ROOT.RooRealVar("bkg_nr_mll", "bkg_nr_mll", 0.01, 0.01, 100)
bkg_mean_mll.setVal(90.9369)
bkg_sigmal_mll.setVal(0.0100032)
bkg_sigmar_mll.setVal(0.709995)
bkg_alphal_mll.setVal(0.709995)
bkg_alphar_mll.setVal(0.505884)
bkg_nl_mll.setVal(0.0533984)
bkg_nr_mll.setVal(0.538716)
bkgrealmll_mll = ROOT.RooCrystalBall("bkgrealmll_mll", "bkgrealmll_mll", mll, mean_mll, sigmal_mll, sigmar_mll, alphal_mll, nl_mll, alphar_mll, nr_mll)

#Background 2d realmll model: bkgrealmll_mll_met_2dpdf = bkgrealmll_met * bkgrealmll_mll
bkgrealmll_mll_met_2dpdf = ROOT.RooProdPdf("bkgrealmll_mll_met_2dpdf", "bkgrealmll_mll_met_2dpdf", [bkgrealmll_mll, bkgrealmll_met])

#Overall 2D bkg model: bkgtot_mll_met_2dpdf = bkgfakemll_mll_met_2dpdf + ratio_realmll * bkgrealmll_mll_met_2dpdf
ratio_realmll = ROOT.RooRealVar("ratio_realmll", "ratio_realmll", 0.44, 0, 1)
ratio_realmll.setVal(0.1)
# TODO: TEMP: IGNORE PRESENCE OF REAL MLL BACKGROUND while we figure out how to do the falling exponential 
bkgtot_mll_met_2dpdf = ROOT.RooAddPdf("bkgtot_mll_met_2dpdf", "bkgtot_mll_met_2dpdf", [bkgrealmll_mll_met_2dpdf, bkgfakemll_mll_met_2dpdf], [ratio_realmll])
# bkgtot_mll_met_2dpdf = bkgfakemll_mll_met_2dpdf

dPdf["bkg"] = {}
for varInfo in variablesInfo:
    dPdf["bkg"][f"hist_{varInfo[0]}"] = {} 
dPdf["bkg"]["hist_mll"]["var"] = mll
dPdf["bkg"]["hist_mll"]["pdf"] = bkgtot_mll_met_2dpdf
dPdf["bkg"]["hist_mll"]["color"] = ROOT.TColor.GetColor("#964a8b") # purple
dPdf["bkg"]["hist_mll"]["label"] = "Background 2D fit: 1D proj"
dPdf["bkg"]["hist_mll"]["name"] = "bkg1Dproj"
dPdf["bkg"]["hist_mll"]["linestyle"] = 1
dPdf["bkg"]["hist_mll"]["ymax"] = 3


dPdf["bkg"]["hist_met"]["var"] = met
dPdf["bkg"]["hist_met"]["pdf"] = bkgtot_mll_met_2dpdf
dPdf["bkg"]["hist_met"]["color"] = ROOT.TColor.GetColor("#964a8b") # purple
dPdf["bkg"]["hist_met"]["label"] = "Background 2D fit: 1D proj"
dPdf["bkg"]["hist_met"]["name"] = "bkg1Dproj"
dPdf["bkg"]["hist_met"]["linestyle"] = 1
dPdf["bkg"]["hist_met"]["ymax"] = 5

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
        leg.SetHeader("Z(ll)H(bb): results of 2D fit")

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

        outdir = f"/eos/user/s/skkwan/www/higgsino/studies/mll-MET-fit-2D"
        sampleName = d[key]["name"]
        canv.SaveAs(f"{sampleName}-{variable}.pdf")
        canv.SaveAs(f"{sampleName}-{variable}.png")
        os.system(f"mv {sampleName}-{variable}.* {outdir}")

        del(canv)
        del(leg)

