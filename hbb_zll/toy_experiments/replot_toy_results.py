import ROOT
import cmsstyle as CMS
import os
import sys
import re
import argparse
from plot_utils import *

if __name__ == "__main__":
    """
    Replot toy experiment results from an existing results file.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("filename")
    parser.add_argument("--ratio-peaking-in", type=float, default=0.088,
                        help="Injected ratio_peaking value (default: 0.088)")
    parser.add_argument("--conv", action='store_true', help="Plot only toy experiments that converged")
    args = parser.parse_args()

    filename = args.filename
    basename = os.path.basename(filename)

    match = re.match(r"toy_results_(\d+)_(\d+)_nExp_(\d+)_nsig_(\d+)_nbkg_(\d+)(\_(\w.*))?.root", basename)
    if not match:
        print(f"Cannot parse filename: {basename}")
        sys.exit(1)

    m1, m2           = int(match.group(1)), int(match.group(2))
    n_experiments    = int(match.group(3))
    n_sig_in         = int(match.group(4))
    n_bkg_in         = int(match.group(5))
    ratio_peaking_in = args.ratio_peaking_in

    n_toys_plotted = 0

    fixed_r = False
    if "fixed_r" in args.filename:
        fixed_r = True

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

        success = (result.status() == 0 and result.covQual() == 3)


        for pname in ["n_sig", "n_bkg", "ratio_peaking"]:
            par = result.floatParsFinal().find(pname)
            if par and par.getError() > 0:
                print("Before checking: args.conv: ", args.conv, "success: ", success)
                if (args.conv and success) or not args.conv:
                    pulls[pname].append((par.getVal() - injected[pname]) / par.getError())
                    fit_vals[pname].append(par.getVal())
                    if pname == "n_sig":
                        n_toys_plotted += 1 


    f.Close()
    n_good = len(pulls["n_sig"])
    print(f"Converged fits (status=0, covQual=3): {n_good} / {n_experiments}")


    eos_dir = "/eos/user/s/skkwan/www/higgsino/studies/mll-MET-fit-2D/toys"

    basename = f"_{m1}_{m2}_nExp_{n_experiments}_nsig_{n_sig_in}_nbkg_{n_bkg_in}"
    if fixed_r:
        basename += "_fixed_r"

    plot_name_suffix = basename 
    if args.conv:
        plot_name_suffix += "_successOnly"

    for pname, pull_vals in pulls.items():
        if not pull_vals:
            print(f"No pulls for {pname}, skipping")
            continue
        make_pull_plot(pull_vals, pname, n_sig_in, n_bkg_in, n_experiments, n_toys_plotted,
                       mass_point=(m1, m2), eos_dir=eos_dir, plot_name_suffix=plot_name_suffix)

    for pname, vals in fit_vals.items():
        if not vals:
            print(f"No fit values for {pname}, skipping")
            continue
        make_distribution_plot(vals, pname, injected[pname], n_sig_in, n_bkg_in, n_experiments, n_toys_plotted,
                               mass_point=(m1, m2), eos_dir=eos_dir, plot_name_suffix=plot_name_suffix)


    dump_data_file = f"/eos/cms/store/group/phys_susy/skkwan/toys/dump_data{basename}.root"
    eos_dir_selected = eos_dir + "/selected_toys"
    print(eos_dir_selected)
    n_sig_min = -1.99*(n_sig_in + n_bkg_in)
    n_sig_max = -1*(n_sig_in + n_bkg_in)
    find_and_plot_selected_toy(filename, dump_data_file,
                                    n_sig_min, n_sig_max,
                                    m1=650, m2=1,
                                    eos_dir=eos_dir_selected, plot_name_suffix=basename, fixed_r=fixed_r)

