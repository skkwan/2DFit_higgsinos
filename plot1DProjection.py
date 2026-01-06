import os, ROOT
import cmsstyle as CMS

ROOT.gInterpreter.Declare("""
float compute_deltaPhi(float phi1, float phi2, const double c = M_PI) {
    auto r = std::fmod(phi2 - phi1, 2.0 * c);
    if (r < -c) {
    r += 2.0 * c;
    }
    else if (r > c) {
    r -= 2.0 * c;
    }
    return r;
}
""")

ROOT.gInterpreter.Declare("""
float compute_deltaR(float eta1, float eta2, float phi1, float phi2) {
    auto deltaR = sqrt(pow(eta1 - eta2, 2) + pow(compute_deltaPhi(phi1, phi2), 2));
    return deltaR;
}
""")

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
    # ["met", "MET / GeV", 40, 0., 400.], 
    # ["weightXyear", "WeightXyear", 40, 0., 10.],
]

ch1, df1 = getTChainRDF(block1, "tree")
ch2, df2 = getTChainRDF(block2, "tree")
# df = df.Define("gen_deltaPhi_ll_ptmiss", "compute_deltaPhi(gen_leps_p4.Phi(), gen_p4_ptmiss.Phi())")

d = {}
dPdf = {}

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

mean_mgg.setVal(124.943)
sigmal_mgg.setVal(1.52331)
sigmar_mgg.setVal(1.18388)
alphal_mgg.setVal(1.66588)
alphar_mgg.setVal(1.17512)
nl_mgg.setVal(2.5573)
alphar_mgg.setVal(1.17512)
nr_mgg.setVal(32.2549)
mgg = ROOT.RooRealVar("mgg", "mgg", 60, 260) 
sig_dcb_mgg = ROOT.RooCrystalBall("sig_dcb_mgg", "sig_dcb_mgg", mgg, mean_mgg, sigmal_mgg, sigmar_mgg, alphal_mgg, nl_mgg, alphar_mgg, nr_mgg)

dPdf["signal"] = {}
dPdf["signal"]["var"] = mgg 
dPdf["signal"]["pdf"] = sig_dcb_mgg
dPdf["signal"]["color"] = ROOT.TColor.GetColor("#f89c20") # light orange
dPdf["signal"]["label"] = "Signal 2D fit projection"
# df2 = df2.DefaultValueFor("GenModel_TChiZH_ZToLL_500_300", False)
# df2 = df2.Filter("GenModel_TChiZH_ZToLL_500_300", "FullSIM: select (500, 300)")
for key in d:
    d[key]["dataframe"].Report().Print()
# report1 = df1.Report()
# report2 = df2.Report()
# report1.Print()
# report2.Print()

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

    # h1 = df1.Histo1D(("var", "var", nBins, xmin, xmax), variable).GetValue()
    # h2 = df2.Histo1D(("var", "var", nBins, xmin, xmax), variable).GetValue()
    # h1 = addOverflow(h1)
    # h2 = addOverflow(h2)
    # normalizeHist(h1)
    # normalizeHist(h2)
    # h1.SetLineColor(ROOT.TColor.GetColor("#ffa90e"))
    # h2.SetLineColor(ROOT.TColor.GetColor("#3f90da"))
    # h1.SetMarkerColor(0)
    # h1.SetLineWidth(2)
    # h2.SetMarkerColor(0)
    # h2.SetLineWidth(2)
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
    leg.SetHeader("bb#gamma#gamma 2D fit: projection to 1D")
    for key in d:
        leg.AddEntry(d[key][f"hist_{variable}"], d[key]["label"], "l")
    # leg.AddEntry(h1, f"(500, 300) GeV", "l")
    # leg.AddEntry(h2, f"(700, 1) GeV", "l")

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
    # CMS.cmsObjectDraw(h1, "HIST")
    # CMS.cmsObjectDraw(h2, "HIST SAME")

    for key in dPdf: 
        plot = dPdf[key]["var"].frame() 
        # dotted line
        # See all draw options: https://root.cern.ch/doc/v636/classRooAbsPdf.html#aa0f2f98d89525302a06a1b7f1b0c2aa6 
        dPdf[key]["pdf"].plotOn(plot, ROOT.RooFit.Name("signal1Dfit"), ROOT.RooFit.LineColor(dPdf[key]["color"]), ROOT.RooFit.LineWidth(2), ROOT.RooFit.LineStyle(2), ROOT.RooFit.MarkerSize(0))
        leg.AddEntry(plot.findObject("signal1Dfit"), dPdf[key]["label"])
        plot.Draw("SAME")

    CMS.cmsObjectDraw(leg)

    CMS.UpdatePad(canv)

    outdir = f"/eos/user/s/skkwan/www/higgsino/studies/from-meraj-fit-2D"
    sampleName = "input"
    canv.SaveAs(f"{sampleName}-{variable}.pdf")
    canv.SaveAs(f"{sampleName}-{variable}.png")
    os.system(f"mv {sampleName}-{variable}.* {outdir}")

    del(canv)
    del(leg)

