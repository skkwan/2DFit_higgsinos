import argparse
import glob
import os
import ROOT

parser = argparse.ArgumentParser(description="Hadd ntuples for 2D fit")
parser.add_argument("-m", "--mode", required=True, choices=["background", "signal"],
                    help="'background' to hadd bkg ntuples, 'signal' to hadd all signal mass points")
args = parser.parse_args()

basedir = "/eos/cms/store/group/phys_susy/skkwan/condorHistogramming/2026-02-25-00h42m-2018-dataMC-with-SR-ntuples"
signal_basedir = "/eos/cms/store/group/phys_susy/skkwan/condorHistogramming/2026-06-06-00h18m-2018-sample-signal-points"

tree_name = "event_tree"

### BACKGROUNDS
samples = {
    "DYJets": [
            "DYJetsToLL_M-50",
            "DYJetsToLL_M-50_HT-70to100",
            "DYJetsToLL_M-50_HT-100to200",
            "DYJetsToLL_M-50_HT-200to400",
            "DYJetsToLL_M-50_HT-400to600",
            "DYJetsToLL_M-50_HT-600to800",
            "DYJetsToLL_M-50_HT-800to1200",
            "DYJetsToLL_M-50_HT-1200to2500",
            "DYJetsToLL_M-50_HT-2500toInf",
        ],
        "WJets": [
            "WJetsToLNu",
        ],
        "ttbar": [
            "TTTo2L2Nu",
            ],
		"TTZ_peak": [
			"TTZToLLNuNu_M-10_peak_mll",
		],
		"TTZ_nonpeak": [
			"TTZToLLNuNu_M-10_nonpeak_mll",
		],
		"TTW": [
			"TTWJetsToLNu",
		],
        "WH": [
            "WminusH_HToBB_WToLNu_M-125",
            "WplusH_HToBB_WToLNu_M-125",
		],
		"ZH": [
            "ZH_HToBB_ZToLL_M-125",
        ],
        "WZ": [
            "WZTo3LNu",
            "WZTo2Q2L",
			"WZTo2Q2Nu",
		],
		"WW": [
            "WWTo2L2Nu",
        ],
        "ZZ": [
            "ZZTo2L2Nu",
            "ZZTo2Q2L",
            "ZZTo2Q2Nu",
			"ZZTo4L",
        ]
}

peaking_samples = ["DYJets", "TTZ_peak"]


def count_yield(class_list):
    """Return total weighted yield for a list of sample classes (e.g. ["DYJets", "TTZ_peak"]).
    mm ntuples are weighted by weight_nominal_mm, ee ntuples by weight_nominal_ee."""
    total = 0.0
    for cls in class_list:
        if cls not in samples:
            raise ValueError(f"Unknown sample class: {cls}")
        for s in samples[cls]:
            mm_files = glob.glob(f"{basedir}/{s}/snapshot*mm_SR_mll_MET_fit_scheme.root")
            ee_files = glob.glob(f"{basedir}/{s}/snapshot*ee_SR_mll_MET_fit_scheme.root")
            if mm_files:
                total += ROOT.RDataFrame(tree_name, mm_files).Sum("weight_nominal_mm").GetValue()
            if ee_files:
                total += ROOT.RDataFrame(tree_name, ee_files).Sum("weight_nominal_ee").GetValue()
    return total


def harmonize_weight(in_file, channel):
    """Define weight_nominal from weight_nominal_{channel} and snapshot to a _fixed.root file."""
    out_file = in_file.replace(".root", "_fixed.root")
    ROOT.RDataFrame(tree_name, in_file) \
        .Define("weight_nominal", f"weight_nominal_{channel}") \
        .Snapshot(tree_name, out_file)
    return out_file


if args.mode == "background":
    # peaking_yield = count_yield(["DYJets", "TTZ_peak"])
    # nonpeaking_yield = count_yield(["WJets", "ttbar", "TTZ_nonpeak", "TTW", "WH", "ZH", "WZ", "WW", "ZZ"])
    # frac = (peaking_yield)/(peaking_yield + nonpeaking_yield)
    # print(f"Counted {peaking_yield} events (weighted) in peaking background and {nonpeaking_yield} (weighted) in non-peaking background. So the peaking background is {frac} of the total background")

    #### BACKGROUND: SEPARATE PEAKING AND NON-PEAKING
    bkg_list_peaking_mm = []
    bkg_list_peaking_ee = []
    bkg_list_nonpeak_mm = []
    bkg_list_nonpeak_ee = []

    for group in samples:
        for s in samples[group]:
            if group in peaking_samples:
                print("Peaking: ", group, s)
                for ntuple in glob.glob(f"{basedir}/{s}/snapshot*mm_SR_mll_MET_fit_scheme.root"):
                    bkg_list_peaking_mm.append(ntuple)
                for ntuple in glob.glob(f"{basedir}/{s}/snapshot*ee_SR_mll_MET_fit_scheme.root"):
                    bkg_list_peaking_ee.append(ntuple)
            else:
                print("Non-peaking: ", group, s)
                for ntuple in glob.glob(f"{basedir}/{s}/snapshot*mm_SR_mll_MET_fit_scheme.root"):
                    bkg_list_nonpeak_mm.append(ntuple)
                for ntuple in glob.glob(f"{basedir}/{s}/snapshot*ee_SR_mll_MET_fit_scheme.root"):
                    bkg_list_nonpeak_ee.append(ntuple)

    # Hadd mm channel PEAKING
    hadd_mm = f"hadd -f -j -k {basedir}/backgrounds_peaking_mm.root"
    for b in bkg_list_peaking_mm:
        hadd_mm += f" {b}"
    print(hadd_mm)
    os.system(hadd_mm)

    # Hadd ee channel PEAKING
    hadd_ee = f"hadd -f -j -k {basedir}/backgrounds_peaking_ee.root"
    for b in bkg_list_peaking_ee:
        hadd_ee += f" {b}"
    print(hadd_ee)
    os.system(hadd_ee)

    # Hadd mm channel NON-PEAKING
    hadd_mm = f"hadd -f -j -k {basedir}/backgrounds_nonpeak_mm.root"
    for b in bkg_list_nonpeak_mm:
        hadd_mm += f" {b}"
    print(hadd_mm)
    os.system(hadd_mm)

    # Hadd ee channel NON-PEAKING
    hadd_ee = f"hadd -f -j -k {basedir}/backgrounds_nonpeak_ee.root"
    for b in bkg_list_nonpeak_ee:
        hadd_ee += f" {b}"
    print(hadd_ee)
    os.system(hadd_ee)

    # Harmonize weight branch name
    harmonize_weight(f"{basedir}/backgrounds_peaking_mm.root", "mm")
    harmonize_weight(f"{basedir}/backgrounds_peaking_ee.root", "ee")

    hadd_combined = f"hadd -f -j -k backgrounds_peaking.root {basedir}/backgrounds_peaking_mm_fixed.root {basedir}/backgrounds_peaking_ee_fixed.root"
    print(hadd_combined)
    os.system(hadd_combined)
    for f in [f"{basedir}/backgrounds_peaking_mm.root", f"{basedir}/backgrounds_peaking_ee.root",
              f"{basedir}/backgrounds_peaking_mm_fixed.root", f"{basedir}/backgrounds_peaking_ee_fixed.root"]:
        os.remove(f)

    # Harmonize weight branch name
    harmonize_weight(f"{basedir}/backgrounds_nonpeak_mm.root", "mm")
    harmonize_weight(f"{basedir}/backgrounds_nonpeak_ee.root", "ee")

    hadd_combined = f"hadd -f -j -k backgrounds_nonpeak.root {basedir}/backgrounds_nonpeak_mm_fixed.root {basedir}/backgrounds_nonpeak_ee_fixed.root"
    print(hadd_combined)
    os.system(hadd_combined)
    for f in [f"{basedir}/backgrounds_nonpeak_mm.root", f"{basedir}/backgrounds_nonpeak_ee.root",
              f"{basedir}/backgrounds_nonpeak_mm_fixed.root", f"{basedir}/backgrounds_nonpeak_ee_fixed.root"]:
        os.remove(f)

    # ### BACKGROUND: TOTAL
    # bkg_list_mm = []
    # bkg_list_ee = []
    # for group in samples:
    #     if "data" in group:
    #         continue
    #     for s in samples[group]:
    #         # print(s)
    #         for ntuple in glob.glob(f"{basedir}/{s}/snapshot*mm_SR_mll_MET_fit_scheme.root"):
    #             bkg_list_mm.append(ntuple)
    #             print(ntuple)
    #         for ntuple in glob.glob(f"{basedir}/{s}/snapshot*ee_SR_mll_MET_fit_scheme.root"):
    #             bkg_list_ee.append(ntuple)

    # # Hadd mm channel
    # hadd_mm = f"hadd -f -j -k {basedir}/backgrounds_mm.root"
    # for b in bkg_list_mm:
    #     hadd_mm += f" {b}"
    # print(hadd_mm)
    # os.system(hadd_mm)

    # # Hadd ee channel
    # hadd_ee = f"hadd -f -j -k {basedir}/backgrounds_ee.root"
    # for b in bkg_list_ee:
    #     hadd_ee += f" {b}"
    # print(hadd_ee)
    # os.system(hadd_ee)

    # rdf_mm = ROOT.RDataFrame(tree_name, f"{basedir}/backgrounds_mm.root")
    # rdf_mm.Define("weight_nominal", "weight_nominal_mm") \
    #       .Snapshot(tree_name, f"{basedir}/backgrounds_mm_fixed.root")

    # rdf_ee = ROOT.RDataFrame(tree_name, f"{basedir}/backgrounds_ee.root")
    # rdf_ee.Define("weight_nominal", "weight_nominal_ee") \
    #       .Snapshot(tree_name, f"{basedir}/backgrounds_ee_fixed.root")

    # # Combine fixed files into a single file
    # hadd_combined = f"hadd -f -j -k backgrounds_for_2D_fit.root {basedir}/backgrounds_mm_fixed.root {basedir}/backgrounds_ee_fixed.root"
    # print(hadd_combined)
    # os.system(hadd_combined)

elif args.mode == "signal":
    #### SIGNAL: hadd all mass points found in signal_basedir
    mass_point_dirs = sorted(glob.glob(f"{signal_basedir}/TChiZH_*"))
    for mp_dir in mass_point_dirs:
        mp_name = os.path.basename(mp_dir)  # e.g. TChiZH_650_1
        mm_files = sorted(glob.glob(f"{mp_dir}/snapshot_{mp_name}_*mm_SR_mll_MET_fit_scheme.root"))
        ee_files = sorted(glob.glob(f"{mp_dir}/snapshot_{mp_name}_*ee_SR_mll_MET_fit_scheme.root"))
        if not mm_files and not ee_files:
            print(f"No snapshot files found in {mp_dir}, skipping.")
            continue

        fixed_files = []
        if mm_files:
            mm_merged = f"{signal_basedir}/snapshot_{mp_name}_mm.root"
            os.system(f"hadd -f -j -k {mm_merged} " + " ".join(mm_files))
            fixed_files.append(harmonize_weight(mm_merged, "mm"))
        if ee_files:
            ee_merged = f"{signal_basedir}/snapshot_{mp_name}_ee.root"
            os.system(f"hadd -f -j -k {ee_merged} " + " ".join(ee_files))
            fixed_files.append(harmonize_weight(ee_merged, "ee"))

        out_file = f"{signal_basedir}/snapshot_{mp_name}_SR_mll_MET_fit_scheme.root"
        os.system(f"hadd -f -j -k {out_file} " + " ".join(fixed_files))
        for f in ([mm_merged] if mm_files else []) + ([ee_merged] if ee_files else []) + fixed_files:
            os.remove(f)
