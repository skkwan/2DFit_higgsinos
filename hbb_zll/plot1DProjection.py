import os, ROOT
import cmsstyle as CMS

doLog = False

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
    ["mll", "m(ll) / GeV", 60, 60., 120., mll],
    ["met", "MET / GeV", 120, 0., 1200., met], 
    # ["weightXyear", "WeightXyear", 40, 0., 10.],
]

resultsfile = ROOT.TFile.Open("fitresult.root", "READ")
workspace = resultsfile.Get("workspace")
signal_results = resultsfile.Get("signal_result")
bkg_results = resultsfile.Get("bkg_result")

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
variables = ROOT.RooArgSet(mll, met, weightXyear)
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


mean_mll = signal_results.floatParsFinal().find("mean_mll")
sigmal_mll = signal_results.floatParsFinal().find("sigmal_mll")
sigmar_mll = signal_results.floatParsFinal().find("sigmar_mll")
alphal_mll = signal_results.floatParsFinal().find("alphal_mll")
alphar_mll = signal_results.floatParsFinal().find("alphar_mll")
nl_mll = signal_results.floatParsFinal().find("nl_mll")
nr_mll = signal_results.floatParsFinal().find("nr_mll")
mypdf = workspace.pdf("sig_roohistpdf_met")
print("mypdf:", mypdf)
myspline = workspace.function("pdf_of_spline")

# Temporary placeholder
sig_dcb_mll = ROOT.RooCrystalBall("sig_dcb_mll", "sig_dcb_mll", mll, mean_mll, sigmal_mll, sigmar_mll, alphal_mll, nl_mll, alphar_mll, nr_mll)

#Signal 2D model: sigtot_mll_met_2dpdf = sig_smoid_met * sig_dcb_mll
sigtot_mll_met_2dpdf = ROOT.RooProdPdf("sigtot_dcb_mll_dcb_met", "sigtot_dcb_mll_dcb_met", [sig_dcb_mll, myspline])

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

#Background real mll model in met dimension
mu_fakemll_met = bkg_results.floatParsFinal().find("mu_fakemll_met")
b_fakemll_met = bkg_results.floatParsFinal().find("b_fakemll_met")
bkgfakemll_met = ROOT.RooGenericPdf("bkgfakemll_met", "bkgfakemll_met", "1/b_fakemll_met * exp(-(@0 - mu_fakemll_met)/b_fakemll_met - exp(-(@0 - mu_fakemll_met)/b_fakemll_met))",
                        ROOT.RooArgList(met, mu_fakemll_met, b_fakemll_met))  
#Background fakemll model in mll dimension
a_fakemll_mll = bkg_results.floatParsFinal().find("a_fakemll_mll")
bkgfakemll_mll = ROOT.RooExponential("bkgfakemll_mll", "bkgfakemll_mll", mll, a_fakemll_mll)
#Background 2d fakemll model: bkgfakemll_mll_met_2dpdf = bkgfakemll_met * bkgfakemll_mll
bkgfakemll_mll_met_2dpdf = ROOT.RooProdPdf("bkgfakemll_mll_met_2dpdf", "bkgfakemll_mll_met_2dpdf", [bkgfakemll_mll, bkgfakemll_met])

#Background realmll model in met dimension
mu_realmll_met = bkg_results.floatParsFinal().find("mu_realmll_met")
b_realmll_met = bkg_results.floatParsFinal().find("b_realmll_met")
bkgrealmll_met = ROOT.RooGenericPdf("bkgrealmll_met", "bkgrealmll_met", "1/b_realmll_met * exp(-(@0 - mu_realmll_met)/b_realmll_met - exp(-(@0 - mu_realmll_met)/b_realmll_met))",
                        ROOT.RooArgList(met, mu_realmll_met, b_realmll_met))  

# Background realmll in mll dimension
bkg_mean_mll = bkg_results.floatParsFinal().find("bkg_mean_mll")
bkg_sigma_mll = bkg_results.floatParsFinal().find("bkg_sigma_mll")
bkgrealmll_mll = ROOT.RooGaussian("bkg_gaus_mll", "bkg_gaus_mll", mll, bkg_mean_mll, bkg_sigma_mll)

#Background 2d realmll model: bkgrealmll_mll_met_2dpdf = bkgrealmll_met * bkgrealmll_mll
bkgrealmll_mll_met_2dpdf = ROOT.RooProdPdf("bkgrealmll_mll_met_2dpdf", "bkgrealmll_mll_met_2dpdf", [bkgrealmll_mll, bkgrealmll_met])

#Overall 2D bkg model: bkgtot_mll_met_2dpdf = bkgfakemll_mll_met_2dpdf + ratio_realmll * bkgrealmll_mll_met_2dpdf
ratio_realmll = bkg_results.floatParsFinal().find("ratio_realmll")

bkgtot_mll_met_2dpdf = ROOT.RooAddPdf("bkgtot_mll_met_2dpdf", "bkgtot_mll_met_2dpdf", [bkgrealmll_mll_met_2dpdf, bkgfakemll_mll_met_2dpdf], [ratio_realmll])


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

        y_min = 0
        y_max = dPdf[key][f"hist_{variable}"]["ymax"]
        if doLog:
            y_min = 1e-10
            y_max = y_max * 1000
        # print(y_maxima)
        # example: https://cms-analysis.docs.cern.ch/guidelines/plotting/examples/?h=pad#stack-plot-with-cmsstyle
        # documentation: https://cmsstyle.readthedocs.io/en/latest/reference/#cmsstyle.cmsstyle.cmsDiCanvas
        canv = CMS.cmsDiCanvas(variable, x_min=xmin, x_max=xmax, y_min=y_min, y_max=y_max, # max(y_maxima) * 1.25,
                            r_min=0,
                            r_max=2,
                            nameXaxis = xlabel,
                            nameYaxis = 'Shape (A.U.)',
                            nameRatio = 'MC/Pred',
                            square = CMS.kSquare,  iPos=0) # yTitOffset=1.6,  # extraSpace=0.05,
        canv.SetRightMargin(0.05)
        # CMS.SetExtraText("Private work (CMS simulation)")
        # CMS.SetCmsTextFont(52)
        # CMS.SetCmsTextSize(0.75*0.76)
        CMS.UpdatePad(canv)

        y_maxima = []

        # Plot the dataset fitted
        canv.cd(1)
        if doLog:
            ROOT.gPad.SetLogy()

        CMS.UpdatePad(canv)
        d[key]["dataset"].plotOn(frame, ROOT.RooFit.Name(d[key]["name"]),
                                    ROOT.RooFit.LineColor(d[key]["color"]),
                                    ROOT.RooFit.LineWidth(d[key]["linewidth"]),
                                    ROOT.RooFit.MarkerColor(d[key]["color"]),
                                    ROOT.RooFit.MarkerSize(1))
        # Hodgepodge of arguments
        d[key][f"hist_{variable}"] = {}
        d[key][f"hist_{variable}"]["var"] = dPdf[key][f"hist_{variable}"]["var"] # copy from dPdf
        d[key][f"hist_{variable}"]["obj"] = d[key]["dataset"].createHistogram("histo", d[key][f"hist_{variable}"]["var"], ROOT.RooFit.Binning(nBins, xmin, xmax))
        leg.AddEntry(frame.findObject(d[key]["name"]), d[key]["label"])

        # Plot the PDF with the fitted parameters
        # See all draw options: https://root.cern.ch/doc/v636/classRooAbsPdf.html#aa0f2f98d89525302a06a1b7f1b0c2aa6 
        dPdf[key][f"hist_{variable}"]["pdf"].plotOn(frame,
                                                    ROOT.RooFit.Name(dPdf[key][f"hist_{variable}"]["name"]),
                                                    ROOT.RooFit.LineColor(dPdf[key][f"hist_{variable}"]["color"]),
                                                    ROOT.RooFit.LineWidth(2),
                                                    ROOT.RooFit.LineStyle(dPdf[key][f"hist_{variable}"]["linestyle"]),
                                                    ROOT.RooFit.MarkerSize(0))
        dPdf[key][f"hist_{variable}"]["RooCurve"] = frame.getCurve(dPdf[key][f"hist_{variable}"]["name"])
        leg.AddEntry(frame.findObject(dPdf[key][f"hist_{variable}"]["name"]), dPdf[key][f"hist_{variable}"]["label"])
        frame.Draw("SAME")

        # Ratio plot
        canv.cd(2)
        data_ratio = d[key][f"hist_{variable}"]["obj"].Clone()
        prediction = d[key][f"hist_{variable}"]["obj"].Clone()
        # At each point in the TH1F, we need to evaluate the pdf value 
        for i in range(1, data_ratio.GetNbinsX() + 1):
            thisXval = prediction.GetXaxis().GetBinCenter(i)
            thisArgSet = ROOT.RooArgSet(mll, met, weightXyear)
            thisArgSet.setRealValue(variable, thisXval)

            curve = dPdf[key][f"hist_{variable}"]["RooCurve"]
            fitY = curve.interpolate(thisXval, tolerance=0.1)
            # fitX = curve.GetPointX(thisPoint)
            # fitY = curve.GetPointY(thisPoint)
            # print(f"Setting {variable} to {thisXval}")
            if "variable" == "met" and "signal" in key:
                # # RooHistPDF block 
                # print("Signal met: special case")
                # fitY = dPdf[key][f"hist_{variable}"]["name"].getVal(thisArgSet)
                # # Spline block
                fitY = spline.getVal(thisArgSet)

            print(f"{variable}: From {thisXval}: at point {thisXval}, {fitY}, compare to data point {thisXval}, {prediction.GetBinContent(i)}")
            prediction.SetBinContent(i, fitY)

        data_ratio.Divide(prediction)
        CMS.cmsObjectDraw(data_ratio, "E", MarkerStyle=ROOT.kFullCircle)
        unitLine = ROOT.TLine(xmin, 1.0, xmax, 1.0)
        unitLine.SetLineColor(ROOT.kBlack)
        unitLine.SetLineWidth(1)
        unitLine.Draw("SAME")
        canv.cd(1)
        CMS.cmsObjectDraw(leg)

        CMS.UpdatePad(canv)

        outdir = f"/eos/user/s/skkwan/www/higgsino/studies/mll-MET-fit-2D"
        sampleName = d[key]["name"]
        plotname = f"{sampleName}-{variable}"
        if doLog:
            plotname = f"{sampleName}-{variable}-log"
        canv.SaveAs(f"{plotname}.pdf")
        canv.SaveAs(f"{plotname}.png")
        os.system(f"mv {plotname}.* {outdir}")

        del(canv)
        del(leg)

