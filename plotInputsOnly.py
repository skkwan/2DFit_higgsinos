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

variablesInfo = [
    ["mgg", "m_{#gamma#gamma} / GeV", 40, 60., 260.],
    ["met", "MET / GeV", 40, 0., 400.], 
    ["weightXyear", "WeightXyear", 40, 0., 10.],
]

ch1, df1 = getTChainRDF(block1, "tree")
ch2, df2 = getTChainRDF(block2, "tree")

d = {}

colors = [ROOT.TColor.GetColor("#832db6"),  # purple
          ROOT.TColor.GetColor("#bd1f01"),  # red
          ROOT.TColor.GetColor("#717581")]  # grey

d["signal"] = {}
d["signal"]["dataframe"] = df1
d["signal"]["color"] = ROOT.TColor.GetColor("#5790fc") # blue 
d["signal"]["label"] = "Signal"
print("Done getting first chain")

d["bkg"] = {}
d["bkg"]["dataframe"] = df2
d["bkg"]["color"] = ROOT.TColor.GetColor("#9c9ca1") # grey
d["bkg"]["label"] = "Background"


for varInfo in variablesInfo:
    variable = varInfo[0]
    xlabel = varInfo[1]
    nBins = varInfo[2]
    xmin = varInfo[3]
    xmax = varInfo[4]

    y_maxima = []

    for key in d:
        d[key][f"hist_{variable}"] = d[key]["dataframe"].Histo1D(("var", "var", nBins, xmin, xmax), variable).GetValue()
        d[key][f"hist_{variable}"] = addOverflow(d[key][f"hist_{variable}"])
        normalizeHist(d[key][f"hist_{variable}"])
        d[key][f"hist_{variable}"].SetLineColor(d[key]["color"]) 
        d[key][f"hist_{variable}"].SetMarkerColor(0)
        d[key][f"hist_{variable}"].SetLineWidth(2)
        y_maxima.append(d[key][f"hist_{variable}"].GetMaximum())

    legYmin = 0.89 - 0.05 * 4
    legYmax = 0.89
    legXmin = 0.5
    legXmax = 0.5 + 0.4
    if ("eta" in variable) or ("phi" in variable):
        legYmin = 0.40 - 0.05 * 4
        legYmax = 0.40
        legXmin = 0.3
        legXmax = 0.4 + 0.4
    leg = CMS.cmsLeg(legXmin, legYmin, legXmax, legYmax, textSize=0.03)
    leg.SetHeader("bb#gamma#gamma: input variables")
    for key in d:
        leg.AddEntry(d[key][f"hist_{variable}"], d[key]["label"], "l")

    CMS.SetLumi("")

    # print(y_maxima)
    canv = CMS.cmsCanvas(variable, xmin, xmax, 0, max(y_maxima) * 1.25,
                        nameXaxis = xlabel,
                        nameYaxis = 'Shape (A.U.)',
                        square = CMS.kSquare, extraSpace=0.05, yTitOffset=1.6, iPos=0)
    canv.SetRightMargin(0.05)
    # CMS.SetExtraText("Private work (CMS simulation)")
    # CMS.SetCmsTextFont(52)
    # CMS.SetCmsTextSize(0.75*0.76)
    CMS.UpdatePad(canv)

    for key in d:
        CMS.cmsObjectDraw(d[key][f"hist_{variable}"], "HIST SAME")
 
    CMS.cmsObjectDraw(leg)

    CMS.UpdatePad(canv)

    outdir = f"/eos/user/s/skkwan/www/higgsino/studies/from-meraj-fit-2D"
    sampleName = "input"
    canv.SaveAs(f"{sampleName}-{variable}.pdf")
    canv.SaveAs(f"{sampleName}-{variable}.png")
    os.system(f"mv {sampleName}-{variable}.* {outdir}")

    del(canv)
    del(leg)