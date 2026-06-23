import ROOT
import cmsstyle as CMS
import os
import sys
import re
import argparse


def make_pull_plot(pull_vals, param_name, n_sig_in, n_bkg_in, n_experiments, mass_point, eos_dir="."):
    m1, m2 = mass_point

    pull_var = ROOT.RooRealVar("pull", "pull", -5, 5)
    pull_ds  = ROOT.RooDataSet("pull_ds", "pull_ds", ROOT.RooArgSet(pull_var))
    for v in pull_vals:
        if -5 < v < 5:
            pull_var.setVal(v)
            pull_ds.add(ROOT.RooArgSet(pull_var))

    pull_frame = pull_var.frame(ROOT.RooFit.Bins(30), ROOT.RooFit.Range(-5, 5), ROOT.RooFit.Title(""))
    pull_ds.plotOn(pull_frame,
                   ROOT.RooFit.Name("pull_data"),
                   ROOT.RooFit.MarkerColor(ROOT.kBlack),
                   ROOT.RooFit.LineColor(ROOT.kBlack),
                   ROOT.RooFit.MarkerStyle(ROOT.kFullCircle),
                   ROOT.RooFit.MarkerSize(0.8))

    mean_fit  = ROOT.RooRealVar("mean_fit",  "mean_fit",  0,   -3, 3)
    sigma_fit = ROOT.RooRealVar("sigma_fit", "sigma_fit", 1, 0.1, 5)
    gauss     = ROOT.RooGaussian("gauss", "gauss", pull_var, mean_fit, sigma_fit)
    gauss.fitTo(pull_ds, ROOT.RooFit.PrintLevel(-1))
    gauss.plotOn(pull_frame,
                 ROOT.RooFit.Name("gaus_fit"),
                 ROOT.RooFit.LineWidth(1),
                 ROOT.RooFit.LineColor(ROOT.TColor.GetColor("#5790fc")))

    CMS.SetExtraText("Private work")
    CMS.SetCmsText("CMS", font=62, size=0.76)
    CMS.SetLumi(250, unit="fb", run="2018")

    y_max = pull_frame.GetMaximum() * 1.4
    canv = CMS.cmsCanvas(f"canv_{param_name}_pull",
                          x_min=-5, x_max=5, y_min=0, y_max=y_max,
                          nameXaxis=f"({param_name} #minus {param_name}^{{in}}) / #sigma_{{fit}}",
                          nameYaxis="Toys",
                          square=True, extraSpace=0.01, iPos=0.)
    canv.SetRightMargin(0.05)
    CMS.UpdatePad(canv)

    pull_frame.Draw("SAME")

    x = 0.15
    leg = ROOT.TLegend(x, 0.68, x+0.80, 0.88)
    leg.SetNColumns(2)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.035)
    leg.AddEntry(ROOT.nullptr, f"Signal ({m1}, {m2}) GeV", "")
    leg.AddEntry(ROOT.nullptr, f"#mu = {mean_fit.getVal():.3f} #pm {mean_fit.getError():.3f}", "")
    leg.AddEntry(ROOT.nullptr, f"N_{{toys}} = {n_experiments}", "")
    leg.AddEntry(ROOT.nullptr, f"#sigma = {sigma_fit.getVal():.3f} #pm {sigma_fit.getError():.3f}", "")
    leg.AddEntry(ROOT.nullptr, f"n_{{sig}}^{{in}} = {n_sig_in},  n_{{bkg}}^{{in}} = {n_bkg_in}", "")
    leg.AddEntry(pull_frame.findObject("gaus_fit"), "Gaussian fit", "l")
    CMS.cmsObjectDraw(leg)
    CMS.UpdatePad(canv)

    plotname = f"pull_{param_name}"

    canv.SaveAs(f"{plotname}.eps")
    os.system(f"gs -q -dBATCH -dNOPAUSE -dSAFER -dEPSCrop -dPDFSETTINGS=/prepress -sDEVICE=pdfwrite "
              f"-dEmbedAllFonts=true -dSubsetFonts=true -sOutputFile={plotname}.pdf {plotname}.eps && rm {plotname}.eps")

    canv.SaveAs(f"{plotname}.png")
    print(f"Created {plotname}.pdf / .png")
    if eos_dir != ".":
        os.system(f"mv {plotname}.pdf {plotname}.png {eos_dir}/")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("filename")
    parser.add_argument("--ratio-peaking-in", type=float, default=0.088,
                        help="Injected ratio_peaking value (default: 0.088)")
    args = parser.parse_args()

    filename = args.filename
    basename = os.path.basename(filename)

    match = re.match(r"toy_results_(\d+)_(\d+)_nExp_(\d+)_nsig_(\d+)_nbkg_(\d+)\.root", basename)
    if not match:
        print(f"Cannot parse filename: {basename}")
        sys.exit(1)

    m1, m2           = int(match.group(1)), int(match.group(2))
    n_experiments    = int(match.group(3))
    n_sig_in         = int(match.group(4))
    n_bkg_in         = int(match.group(5))
    ratio_peaking_in = args.ratio_peaking_in

    print(f"Mass point: ({m1}, {m2}),  n_experiments={n_experiments},  n_sig_in={n_sig_in},  n_bkg_in={n_bkg_in},  ratio_peaking_in={ratio_peaking_in}")

    f = ROOT.TFile.Open(filename, "READ")

    injected = {"n_sig": n_sig_in, "n_bkg": n_bkg_in, "ratio_peaking": ratio_peaking_in}
    pulls    = {"n_sig": [], "n_bkg": [], "ratio_peaking": []}

    for i in range(n_experiments):
        result = f.Get(f"fitResult_{i}")
        if not result:
            print(f"Warning: fitResult_{i} not found, skipping")
            continue
        for pname in ["n_sig", "n_bkg", "ratio_peaking"]:
            par = result.floatParsFinal().find(pname)
            if par and par.getError() > 0:
                pulls[pname].append((par.getVal() - injected[pname]) / par.getError())

    f.Close()

    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

    eos_dir = "/eos/user/s/skkwan/www/higgsino/studies/mll-MET-fit-2D/toys"

    for pname, pull_vals in pulls.items():
        if not pull_vals:
            print(f"No pulls for {pname}, skipping")
            continue
        make_pull_plot(pull_vals, pname, n_sig_in, n_bkg_in, n_experiments,
                       mass_point=(m1, m2), eos_dir=eos_dir)
