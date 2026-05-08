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

block = [ "backgrounds_CRZ_Zpeak_2018_mm.root" ]
ch, df = getTChainRDF(block, "event_tree")


mll = ROOT.RooRealVar("m_ll", "m_ll", 60, 120) 
weightXyear = ROOT.RooRealVar("weight_nominal_mm", "weight_nominal_mm", -1, 1)

variablesInfo = [
    ["mll", "m(ll) / GeV", 60, 60., 120., mll],
    # ["weightXyear", "WeightXyear", 40, 0., 10.],
]

### Get Z-peak initial fit results
zpeak_file = ROOT.TFile.Open("initial_zPeak_fit_result.root", "READ")
zpeak_results = zpeak_file.Get("zPeak_CRZ_fit_result")
zpeak_mean = zpeak_results.floatParsFinal().find("bkg_Zpeak_mean_mll")
zpeak_sigma = zpeak_results.floatParsFinal().find("bkg_Zpeak_sigma_mll")

# Background realmll in mll dimension
# bkgrealmll_mll = ROOT.RooGaussian("bkg_gaus_mll", "bkg_gaus_mll", mll, zpeak_mean, zpeak_sigma)


peak_mean_mll = zpeak_results.floatParsFinal().find("peak_mean_mll")
peak_sigmal_mll = zpeak_results.floatParsFinal().find("peak_sigmal_mll")
peak_sigmar_mll = zpeak_results.floatParsFinal().find("peak_sigmar_mll")
peak_alphal_mll = zpeak_results.floatParsFinal().find("peak_alphal_mll")
peak_nl_mll = zpeak_results.floatParsFinal().find("peak_nl_mll")
peak_alphar_mll = zpeak_results.floatParsFinal().find("peak_alphar_mll")
peak_nr_mll = zpeak_results.floatParsFinal().find("peak_nr_mll")
peak_dcb_mll = ROOT.RooCrystalBall("peak_dcb_mll", "peak_dcb_mll", mll, peak_mean_mll, peak_sigmal_mll, peak_sigmar_mll, peak_alphal_mll, peak_nl_mll, peak_alphar_mll, peak_nr_mll)


###### Retrive cr data root file ########
crfilepath = 'backgrounds_CRZ_Zpeak_2018_mm.root'
crfile = ROOT.TFile.Open(crfilepath, "READ")
crtree = crfile.Get("event_tree")
variables = ROOT.RooArgSet(mll, weightXyear)
bkgdataset = ROOT.RooDataSet("bkgdataset", "bkgdataset", variables, ROOT.RooFit.Import(crtree), ROOT.RooFit.WeightVar(weightXyear))


d = {}
dPdf = {}

d["bkg"] = {}
d["bkg"]["dataframe"] = df
d["bkg"]["color"] = ROOT.TColor.GetColor("#9c9ca1") # grey
d["bkg"]["label"] = "MC: DYJets + TTZ peaking"
d["bkg"]["linewidth"] = 2
d["bkg"]["dataset"] = bkgdataset 
d["bkg"]["name"] = "background"


dPdf["bkg"] = {}
for varInfo in variablesInfo:
    dPdf["bkg"][f"hist_{varInfo[0]}"] = {} 
dPdf["bkg"]["hist_mll"]["var"] = mll
dPdf["bkg"]["hist_mll"]["pdf"] = peak_dcb_mll
dPdf["bkg"]["hist_mll"]["color"] = ROOT.TColor.GetColor("#e42536") # red
dPdf["bkg"]["hist_mll"]["label"] = "Z-peak fit (DCB)"
dPdf["bkg"]["hist_mll"]["name"] = "bkg_zPeak"
dPdf["bkg"]["hist_mll"]["linestyle"] = 1
dPdf["bkg"]["hist_mll"]["ymax"] = 30000

if doLog:
    dPdf["bkg"]["hist_mll"]["ymax"] = 30000000


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
        leg.SetHeader("Z(ll)H(bb): initial Z-peak background fit")

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
                            nameYaxis = 'Yield',
                            nameRatio = 'MC/Pred',
                            square = CMS.kSquare, iPos=0, extraSpace=0.1) # yTitOffset=1.6,  # extraSpace=0.05,
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
            thisArgSet = ROOT.RooArgSet(mll, weightXyear)
            thisArgSet.setRealValue(variable, thisXval)

            mll.setVal(thisXval)
            bin_width = 1
            fitY = dPdf[key][f"hist_{variable}"]["pdf"].getVal(ROOT.RooArgSet(mll)) * d[key]["dataset"].sumEntries() * 1


            # curve = dPdf[key][f"hist_{variable}"]["RooCurve"]
            # fitY = curve.interpolate(thisXval, tolerance=0.1)



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
        plotname = f"zPeak-{sampleName}-{variable}"
        if doLog:
            plotname = f"zPeak-{sampleName}-{variable}-log"
        canv.SaveAs(f"{plotname}.pdf")
        canv.SaveAs(f"{plotname}.png")
        os.system(f"mv {plotname}.* {outdir}")

        del(canv)
        del(leg)

