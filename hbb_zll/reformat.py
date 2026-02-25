import glob
import os
import ROOT 


basedir = "/eos/cms/store/group/phys_susy/skkwan/condorHistogramming/2026-02-25-00h42m-2018-dataMC-with-SR-ntuples"

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
        "data_obs": [
            "DoubleMuon_Run2018A",
            "DoubleMuon_Run2018B",
            "DoubleMuon_Run2018C",
            "DoubleMuon_Run2018D",
            "EGamma_Run2018A",
            "EGamma_Run2018B",
            "EGamma_Run2018C",
            "EGamma_Run2018D",
			# MuonEG is overwritten in histConsolidation.py
			# "MuonEG_Run2018A",
			# "MuonEG_Run2018B",
        	# "MuonEG_Run2018C",
        	# "MuonEG_Run2018D",
        ],  
        "TTZ": [
            "TTZToLLNuNu_M-10",
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
			# "ZZ",
            "ZZTo2L2Nu",
            "ZZTo2Q2L",
            "ZZTo2Q2Nu",
			"ZZTo4L", 
        ]
}

bkg_list = []

for group in samples:
    if "data" in group:
        continue
    # print(samples[group])
    for s in samples[group]:
        # print(f"{basedir}/{s}/snapshot*.root")
        ntuples = glob.glob(f"{basedir}/{s}/snapshot*mm_SR_mll_MET_fit_scheme.root")
        for ntuple in ntuples:
            # print(ntuple)
            bkg_list.append(ntuple)

# print(bkg_list)

# Build hadd command
hadd_cmd = "hadd -f -j -k backgrounds.root"
for b in bkg_list:
    hadd_cmd += f" {b}"

print(hadd_cmd)
print(">>>> Command not executed, copy and paste and run the above command")