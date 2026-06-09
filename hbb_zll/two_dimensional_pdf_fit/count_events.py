import os
import ROOT

# os.system("hadd -f -j -k signal_plus_background.root ../pdf_fit/snapshot_TChiZH_650_1_SR_mll_MET_fit_scheme.root ../pdf_fit/backgrounds_for_2D_fit.root")

sig_file = ROOT.TFile.Open("../pdf_fit/snapshot_TChiZH_650_1_SR_mll_MET_fit_scheme.root", "READ")
bkg_file = ROOT.TFile.Open("../pdf_fit/backgrounds_for_2D_fit.root", "READ")

def sum_weights(tree):
    return sum(getattr(ev, "weight_nominal") for ev in tree)

sig_tree = sig_file.Get("event_tree")
bkg_tree = bkg_file.Get("event_tree")

n_sig    = sig_tree.GetEntries()
n_bkg    = bkg_tree.GetEntries()
yield_sig = sum_weights(sig_tree)
yield_bkg = sum_weights(bkg_tree)

print(f"TChiZH signal events:    {n_sig}")
print(f"TChiZH signal yield:     {yield_sig:.4f}")
print(f"Background events:       {n_bkg}")
print(f"Background yield:        {yield_bkg:.4f}")
print(f"Total events:            {n_sig + n_bkg}")
print(f"Total yield:             {yield_sig + yield_bkg:.4f}")

sig_file.Close()
bkg_file.Close()