import glob
import os
import ROOT 


basedir = "/eos/cms/store/group/phys_susy/skkwan/condorHistogramming/2026-02-25-00h42m-2018-dataMC-with-SR-ntuples"

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

samples_signal = {
    "snapshot_TChiZH_650_1_SR_mll_MET_fit_scheme": [
        "/eos/cms/store/group/phys_susy/skkwan/condorHistogramming/2026-02-25-00h42m-2018-dataMC-with-SR-ntuples/TChiZH_650_1/snapshot_TChiZH_650_1_cat_0_batch_0_channel_mm_SR_mll_MET_fit_scheme.root",
        "/eos/cms/store/group/phys_susy/skkwan/condorHistogramming/2026-02-25-00h42m-2018-dataMC-with-SR-ntuples/TChiZH_650_1/snapshot_TChiZH_650_1_cat_0_batch_0_channel_ee_SR_mll_MET_fit_scheme.root",
    ]
}

bkg_list_mm = []
bkg_list_ee = []


for group in samples:
    if "data" in group:
        continue
    for s in samples[group]:
        # print(s)
        for ntuple in glob.glob(f"{basedir}/{s}/snapshot*mm_SR_mll_MET_fit_scheme.root"):
            bkg_list_mm.append(ntuple)
            print(ntuple)
        for ntuple in glob.glob(f"{basedir}/{s}/snapshot*ee_SR_mll_MET_fit_scheme.root"):
            bkg_list_ee.append(ntuple)


# Hadd mm channel
hadd_mm = f"hadd -f -j -k {basedir}/backgrounds_mm.root"
for b in bkg_list_mm:
    hadd_mm += f" {b}"
print(hadd_mm)
os.system(hadd_mm)

# Hadd ee channel
hadd_ee = f"hadd -f -j -k {basedir}/backgrounds_ee.root"
for b in bkg_list_ee:
    hadd_ee += f" {b}"
print(hadd_ee)
os.system(hadd_ee)


rdf_mm = ROOT.RDataFrame(tree_name, f"{basedir}/backgrounds_mm.root")
rdf_mm.Define("weight_nominal", "weight_nominal_mm") \
      .Snapshot(tree_name, f"{basedir}/backgrounds_mm_fixed.root")

rdf_ee = ROOT.RDataFrame(tree_name, f"{basedir}/backgrounds_ee.root")
rdf_ee.Define("weight_nominal", "weight_nominal_ee") \
      .Snapshot(tree_name, f"{basedir}/backgrounds_ee_fixed.root")

# Combine fixed files into a single file
hadd_combined = f"hadd -f -j -k backgrounds_for_2D_fit.root {basedir}/backgrounds_mm_fixed.root {basedir}/backgrounds_ee_fixed.root"
print(hadd_combined)
os.system(hadd_combined)

##### SIGNAL
for out_name, files in samples_signal.items():
    print(out_name, files)
    sig_mm_file = next((f for f in files if "_mm_" in f))
    sig_ee_file = next((f for f in files if "_ee_" in f))
    print(sig_mm_file, sig_ee_file)

    # This is all happening in the EOS area
    rdf_sig_mm = ROOT.RDataFrame(tree_name, sig_mm_file)
    rdf_sig_mm.Define("weight_nominal", "weight_nominal_mm") \
              .Snapshot(tree_name, sig_mm_file.replace(".root", "_fixed.root"))

    rdf_sig_ee = ROOT.RDataFrame(tree_name, sig_ee_file)
    rdf_sig_ee.Define("weight_nominal", "weight_nominal_ee") \
              .Snapshot(tree_name, sig_ee_file.replace(".root", "_fixed.root"))

    hadd_sig = f"hadd -f -j -k {out_name}.root"
    hadd_sig += f" {sig_mm_file.replace('.root', '_fixed.root')}"
    hadd_sig += f" {sig_ee_file.replace('.root', '_fixed.root')}"
    print(hadd_sig)
    os.system(hadd_sig)
