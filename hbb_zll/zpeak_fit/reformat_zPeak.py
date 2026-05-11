import glob
import os
import ROOT


basedir = "/eos/cms/store/group/phys_susy/skkwan/condorHistogramming/2026-04-09-01h05m-2018-CRZ-MC-only-for-fit"

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
		"TTZ_peak": [
			"TTZToLLNuNu_M-10_peak_mll",
		],
}

bkg_list_mm = []
bkg_list_ee = []

for group in samples:
    if "data" in group:
        continue
    for s in samples[group]:
        for ntuple in glob.glob(f"{basedir}/{s}/snapshot*_mm_CRZ.root"):
            bkg_list_mm.append(ntuple)
        for ntuple in glob.glob(f"{basedir}/{s}/snapshot*_ee_CRZ.root"):
            bkg_list_ee.append(ntuple)

# Hadd mm channel
hadd_mm = "hadd -f -j -k backgrounds_CRZ_Zpeak_2018_mm.root"
for b in bkg_list_mm:
    hadd_mm += f" {b}"
print(hadd_mm)
os.system(hadd_mm)

# Hadd ee channel
hadd_ee = "hadd -f -j -k backgrounds_CRZ_Zpeak_2018_ee.root"
for b in bkg_list_ee:
    hadd_ee += f" {b}"
print(hadd_ee)
os.system(hadd_ee)

# Create a default branch weight_nominal 
tree_name = "event_tree"

rdf_mm = ROOT.RDataFrame(tree_name, "backgrounds_CRZ_Zpeak_2018_mm.root")
rdf_mm.Define("weight_nominal", "weight_nominal_mm") \
      .Snapshot(tree_name, "backgrounds_CRZ_Zpeak_2018_mm_fixed.root")

rdf_ee = ROOT.RDataFrame(tree_name, "backgrounds_CRZ_Zpeak_2018_ee.root")
rdf_ee.Define("weight_nominal", "weight_nominal_ee") \
      .Snapshot(tree_name, "backgrounds_CRZ_Zpeak_2018_ee_fixed.root")

# Combine fixed files into a single file
hadd_combined = "hadd -f -j -k backgrounds_CRZ_Zpeak_2018.root backgrounds_CRZ_Zpeak_2018_mm_fixed.root backgrounds_CRZ_Zpeak_2018_ee_fixed.root"
print(hadd_combined)
os.system(hadd_combined)

# Clean up
os.system("rm *mm*.root *ee*.root")