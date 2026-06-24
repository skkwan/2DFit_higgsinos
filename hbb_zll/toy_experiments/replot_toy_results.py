import ROOT
import cmsstyle as CMS
import os
import sys
import re
import argparse

def addOverflow(h: ROOT.TH1F) -> ROOT.TH1F:
    """
    Add overflow to a histogram
    """
    h.SetBinContent(h.GetNbinsX(), h.GetBinContent(h.GetNbinsX()) + h.GetBinContent(h.GetNbinsX() + 1))
    return h

def make_pull_plot(pull_vals, param_name, n_sig_in, n_bkg_in, n_experiments, mass_point, eos_dir="."):
    """
    Make pull plots
    """
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
                 ROOT.RooFit.LineWidth(2),
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


def make_distribution_plot(fit_vals, param_name, injected_val, n_sig_in, n_bkg_in, n_experiments, mass_point, eos_dir="."):
    """
    Make distribution plots of the fitted values in the toy experiments
    """
    m1, m2 = mass_point

    mean_data = sum(fit_vals) / len(fit_vals)
    std_data  = (sum((v - mean_data) ** 2 for v in fit_vals) / len(fit_vals)) ** 0.5
    spread    = max(std_data, abs(mean_data) * 0.01, 1e-9)
    x_min     = mean_data - 4.5 * spread
    x_max     = mean_data + 4.5 * spread
    # Override default distribution
    if param_name == "ratio_peaking":
        x_min, x_max = -0.2, 1.0

    val_var = ROOT.RooRealVar("val", param_name, x_min, x_max)
    val_ds  = ROOT.RooDataSet("val_ds", "val_ds", ROOT.RooArgSet(val_var))
    for v in fit_vals:
        if x_min < v < x_max:
            val_var.setVal(v)
            val_ds.add(ROOT.RooArgSet(val_var))

    val_frame = val_var.frame(ROOT.RooFit.Bins(30), ROOT.RooFit.Range(x_min, x_max), ROOT.RooFit.Title(""))
    val_ds.plotOn(val_frame,
                  ROOT.RooFit.Name("val_data"),
                  ROOT.RooFit.MarkerColor(ROOT.kBlack),
                  ROOT.RooFit.LineColor(ROOT.kBlack),
                  ROOT.RooFit.MarkerStyle(ROOT.kFullCircle),
                  ROOT.RooFit.MarkerSize(0.8))

    mean_fit  = ROOT.RooRealVar("mean_fit",  "mean_fit",  mean_data, x_min, x_max)
    sigma_fit = ROOT.RooRealVar("sigma_fit", "sigma_fit", spread, spread * 0.01, spread * 20)
    gauss     = ROOT.RooGaussian("gauss", "gauss", val_var, mean_fit, sigma_fit)
    gauss.fitTo(val_ds, ROOT.RooFit.PrintLevel(-1))
    gauss.plotOn(val_frame,
                 ROOT.RooFit.Name("gaus_fit"),
                 ROOT.RooFit.LineWidth(2),
                 ROOT.RooFit.LineColor(ROOT.TColor.GetColor("#92dadd")))

    CMS.SetExtraText("Private work")
    CMS.SetCmsText("CMS", font=62, size=0.76)
    CMS.SetLumi(250, unit="fb", run="2018")

    display = {"n_sig": "Best-fit n_{sig}", "n_bkg": "Best-fit n_{bkg}", "ratio peaking": "Best-fit ratio"}
    label   = display.get(param_name, param_name)

    y_max = val_frame.GetMaximum() * 1.4
    canv = CMS.cmsCanvas(f"canv_{param_name}_dist",
                          x_min=x_min, x_max=x_max, y_min=0, y_max=y_max,
                          nameXaxis=label,
                          nameYaxis="Toys",
                          square=True, extraSpace=0.01, iPos=0.)
    canv.SetRightMargin(0.05)
    CMS.UpdatePad(canv)

    val_frame.Draw("SAME")

    x = 0.15
    leg = ROOT.TLegend(x, 0.68, x + 0.80, 0.88)
    leg.SetNColumns(2)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.035)
    leg.AddEntry(ROOT.nullptr, f"Signal ({m1}, {m2}) GeV", "")
    leg.AddEntry(ROOT.nullptr, f"#mu = {mean_fit.getVal():.3f} #pm {mean_fit.getError():.3f}", "")
    leg.AddEntry(ROOT.nullptr, f"N_{{toys}} = {n_experiments}", "")
    leg.AddEntry(ROOT.nullptr, f"#sigma = {sigma_fit.getVal():.3f} #pm {sigma_fit.getError():.3f}", "")
    leg.AddEntry(ROOT.nullptr, f"{label}^{{in}} = {injected_val}", "")
    leg.AddEntry(val_frame.findObject("gaus_fit"), "Gaussian fit", "l")
    CMS.cmsObjectDraw(leg)
    CMS.UpdatePad(canv)

    plotname = f"dist_{param_name}"

    canv.SaveAs(f"{plotname}.eps")
    os.system(f"gs -q -dBATCH -dNOPAUSE -dSAFER -dEPSCrop -dPDFSETTINGS=/prepress -sDEVICE=pdfwrite "
              f"-dEmbedAllFonts=true -dSubsetFonts=true -sOutputFile={plotname}.pdf {plotname}.eps && rm {plotname}.eps")

    canv.SaveAs(f"{plotname}.png")
    print(f"Created {plotname}.pdf / .png")
    if eos_dir != ".":
        os.system(f"mv {plotname}.pdf {plotname}.png {eos_dir}/")


def find_and_plot_selected_toy(results_file, dump_data_file,
                               n_sig_min, n_sig_max,
                               m1=650, m2=1,
                               eos_dir="."):
    """
    Find the first toy experiment where n_sig is in [n_sig_min, n_sig_max],
    reconstruct the model, and plot the fit with '_selected' in the output name.
    """
    import sys as _sys
    _sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import toy as _toy

    # Scan fit results to find the first matching experiment
    f_results = ROOT.TFile.Open(results_file, "READ")
    n_total = f_results.GetListOfKeys().GetEntries()
    selected_i = None
    for i in range(n_total):
        result = f_results.Get(f"fitResult_{i}")
        if not result:
            continue
        par = result.floatParsFinal().find("n_sig")
        if par and n_sig_min <= par.getVal() <= n_sig_max:
            selected_i = i
            break

    if selected_i is None:
        print(f"No toy found with n_sig in [{n_sig_min}, {n_sig_max}]")
        f_results.Close()
        return

    selected_result = f_results.Get(f"fitResult_{selected_i}")
    print(f"Found toy {selected_i}: n_sig = {selected_result.floatParsFinal().find('n_sig').getVal():.2f}")

    # Load the generated dataset (keep file open until after plotting)
    f_data = ROOT.TFile.Open(dump_data_file, "READ")
    gen_data = f_data.Get(f"genData_{selected_i}")
    if not gen_data:
        print(f"genData_{selected_i} not found in {dump_data_file}")
        f_results.Close()
        f_data.Close()
        return

    # Reconstruct the model.
    # get_signal_model / get_background_model in toy.py reference mll as a module-level
    # global, so inject it here before calling them.
    mll = ROOT.RooRealVar("m_ll", "m_ll", 60, 120)
    _toy.mll = mll

    sig_model, met, sig_components = _toy.get_signal_model(m1, m2)
    bkg_model, bkg_components, ratio_peaking = _toy.get_background_model(met)

    n_sig = ROOT.RooRealVar("n_sig", "n_sig",  0, -110, 110)
    n_bkg = ROOT.RooRealVar("n_bkg", "n_bkg", 22, -110, 110)
    model = ROOT.RooAddPdf("total_pdf", "total_pdf",
                           ROOT.RooArgList(sig_model, bkg_model),
                           ROOT.RooArgList(n_sig, n_bkg))

    # Sync floating parameters to fit-result values so component fractions are correct
    for pname, var in [("n_sig", n_sig), ("n_bkg", n_bkg), ("ratio_peaking", ratio_peaking)]:
        par = selected_result.floatParsFinal().find(pname)
        if par:
            var.setVal(par.getVal())

    _toy.plot_toy_fit(gen_data, model, selected_result, met, mll,
                     mass_point=(m1, m2),
                     plotname_prefix=f"toy_fit_i{selected_i}_selected",
                     eos_dir=eos_dir)

    f_data.Close()
    f_results.Close()


if __name__ == "__main__":
    """
    Replot toy experiment results from an existing results file.
    """
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

    injected  = {"n_sig": n_sig_in, "n_bkg": n_bkg_in, "ratio_peaking": ratio_peaking_in}
    pulls     = {"n_sig": [], "n_bkg": [], "ratio_peaking": []}
    fit_vals  = {"n_sig": [], "n_bkg": [], "ratio_peaking": []}

    for i in range(n_experiments):
        result = f.Get(f"fitResult_{i}")
        if not result:
            print(f"Warning: fitResult_{i} not found, skipping")
            continue
        for pname in ["n_sig", "n_bkg", "ratio_peaking"]:
            par = result.floatParsFinal().find(pname)
            if par and par.getError() > 0:
                pulls[pname].append((par.getVal() - injected[pname]) / par.getError())
                fit_vals[pname].append(par.getVal())

    f.Close()

    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

    eos_dir = "/eos/user/s/skkwan/www/higgsino/studies/mll-MET-fit-2D/toys"

    for pname, pull_vals in pulls.items():
        if not pull_vals:
            print(f"No pulls for {pname}, skipping")
            continue
        make_pull_plot(pull_vals, pname, n_sig_in, n_bkg_in, n_experiments,
                       mass_point=(m1, m2), eos_dir=eos_dir)

    for pname, vals in fit_vals.items():
        if not vals:
            print(f"No fit values for {pname}, skipping")
            continue
        make_distribution_plot(vals, pname, injected[pname], n_sig_in, n_bkg_in, n_experiments,
                               mass_point=(m1, m2), eos_dir=eos_dir)


    dump_data_file = f"/eos/cms/store/group/phys_susy/skkwan/toys/dump_data_{m1}_{m2}_nExp_{n_experiments}_nsig_{n_sig_in}_nbkg_{n_bkg_in}.root"
    n_sig_min = -30
    n_sig_max = -20
    find_and_plot_selected_toy(filename, dump_data_file,
                                    n_sig_min, n_sig_max,
                                    m1=650, m2=1,
                                    eos_dir=eos_dir)