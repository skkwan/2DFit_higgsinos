#!/usr/bin/env python3
"""
toy_study.py — 2D (m_ll, MET) toy MC study for the ZH higgsino search.

Uses RooMCStudy.generateAndFit to generate and fit toy datasets from the same
fixed best-fit PDFs as the 2D analysis:
  - Signal    mll: double-sided Crystal Ball (from individual signal fit result)
  - Signal    MET: spline PDF from workspace  (no free parameters)
  - Bkg non-peak  mll x MET: Exponential x Gumbel  (from background fit result)
  - Bkg peaking   mll x MET: Crystal Ball x Gumbel  (from zpeak / bkg fit results)

True yields and r are fixed:
  n_bkg = 19.24,  r = 0.088  (peaking fraction, constant in generation and fit)
  n_sig = 2.32  in --mode sb  (signal + background, mass point 650/1 GeV)
  n_sig = 0.0   in --mode b   (background-only)

Generated datasets are retained (keepGenData=True).  Two output canvases:
  1. toy_study_{mode}_{n_toys}toys.{png,pdf}
       Fitted n_sig / n_bkg distributions (with Gaussian overlay) and their pulls.
  2. toy_study_overlay_{mode}_{n_toys}toys.{png,pdf}
       Box-and-whisker plot of per-bin event counts across all generated toys.

Usage:
  python3 toy_study.py [--n-toys N] [--mode {sb,b}] [--seed S]
"""

import argparse
import math
import ROOT
from ROOT import RooFit as RF
import os

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

# ─────────────────────────────────────────────────────────────────────────────
# Args
# ─────────────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description="2D toy MC study for ZH higgsino search")
parser.add_argument("--n-toys", type=int, default=100, help="Number of toy experiments")
parser.add_argument("--mode",   choices=["sb", "b"], default="sb",
                    help="sb = signal+background (n_sig=2.32), b = background-only (n_sig=0)")
parser.add_argument("--seed",   type=int, default=42,  help="Random seed")
args = parser.parse_args()

print(f"Running {args.n_toys} toys in '{args.mode}' mode  (mass point 650/1 GeV, seed={args.seed})")

# ─────────────────────────────────────────────────────────────────────────────
# Fixed true yields
# ─────────────────────────────────────────────────────────────────────────────
TRUE_N_BKG = 19.24
TRUE_R     = 0.088   # ratio_peaking: peaking fraction of total background
TRUE_N_SIG = 2.32 if args.mode == "sb" else 0.0

# ─────────────────────────────────────────────────────────────────────────────
# Load fit result files
# ─────────────────────────────────────────────────────────────────────────────
base = os.path.dirname(os.path.abspath(__file__))

sig_file   = ROOT.TFile.Open(os.path.join(base,
    "../individual_pdf_fits/individual_fit_results/fitresult_signal_650_1.root"))
bkg_file   = ROOT.TFile.Open(os.path.join(base,
    "../individual_pdf_fits/individual_fit_results/fitresult_background_all_except_ZPeak.root"))
zpeak_file = ROOT.TFile.Open(os.path.join(base,
    "../zpeak_fit/initial_zPeak_fit_result.root"))

for label, f in [("signal", sig_file), ("background", bkg_file), ("zpeak", zpeak_file)]:
    if not f or f.IsZombie():
        raise FileNotFoundError(f"Could not open {label} fit result file")

workspace         = sig_file.Get("workspace_650_1")
sig_result        = sig_file.Get("sig_result")
bkg_np_met_result = bkg_file.Get("bkg_nonpeak_result_met")
bkg_np_mll_result = bkg_file.Get("bkg_nonpeak_result_mll")
bkg_pk_met_result = bkg_file.Get("bkg_peaking_result_met")
zpeak_result      = zpeak_file.Get("zPeak_CRZ_fit_result")

# ─────────────────────────────────────────────────────────────────────────────
# Observables
# Use workspace's met — pdf_of_spline is built against it.
# ─────────────────────────────────────────────────────────────────────────────
met = workspace.var("met")
mll = ROOT.RooRealVar("m_ll", "m_ll", 60, 120)
obs_set = ROOT.RooArgSet(mll, met)

# ─────────────────────────────────────────────────────────────────────────────
# Shape parameters — all fixed at best-fit values
# ─────────────────────────────────────────────────────────────────────────────
_sig_pars   = sig_result.floatParsFinal()
_bgnp_met_p = bkg_np_met_result.floatParsFinal()
_bgnp_mll_p = bkg_np_mll_result.floatParsFinal()
_bgpk_met_p = bkg_pk_met_result.floatParsFinal()
_zpeak_pars = zpeak_result.floatParsFinal()

mu_nonpeak_met = _bgnp_met_p.find("mu_nonpeak_met")
b_nonpeak_met  = _bgnp_met_p.find("b_nonpeak_met")
a_nonpeak_mll  = _bgnp_mll_p.find("a_nonpeak_mll")
mu_peaking_met = _bgpk_met_p.find("mu_peaking_met")
b_peaking_met  = _bgpk_met_p.find("b_peaking_met")
zp_mean_mll    = _zpeak_pars.find("peak_mean_mll")
zp_sigmal_mll  = _zpeak_pars.find("peak_sigmal_mll")
zp_sigmar_mll  = _zpeak_pars.find("peak_sigmar_mll")
zp_alphal_mll  = _zpeak_pars.find("peak_alphal_mll")
zp_alphar_mll  = _zpeak_pars.find("peak_alphar_mll")
zp_nl_mll      = _zpeak_pars.find("peak_nl_mll")
zp_nr_mll      = _zpeak_pars.find("peak_nr_mll")
sig_mean_mll   = _sig_pars.find("mean_mll_650_1")
sig_sigmal_mll = _sig_pars.find("sigmal_mll_650_1")
sig_sigmar_mll = _sig_pars.find("sigmar_mll_650_1")
sig_alphal_mll = _sig_pars.find("alphal_mll_650_1")
sig_alphar_mll = _sig_pars.find("alphar_mll_650_1")
sig_nl_mll     = _sig_pars.find("nl_mll_650_1")
sig_nr_mll     = _sig_pars.find("nr_mll_650_1")

for v in [mu_nonpeak_met, b_nonpeak_met, a_nonpeak_mll,
          mu_peaking_met, b_peaking_met,
          zp_mean_mll, zp_sigmal_mll, zp_sigmar_mll,
          zp_alphal_mll, zp_alphar_mll, zp_nl_mll, zp_nr_mll,
          sig_mean_mll, sig_sigmal_mll, sig_sigmar_mll,
          sig_alphal_mll, sig_alphar_mll, sig_nl_mll, sig_nr_mll]:
    v.setConstant(True)

print("Shape parameters (fixed):")
for v in [mu_nonpeak_met, b_nonpeak_met, a_nonpeak_mll,
          mu_peaking_met, b_peaking_met,
          zp_mean_mll, zp_sigmal_mll, zp_sigmar_mll,
          zp_alphal_mll, zp_alphar_mll, zp_nl_mll, zp_nr_mll,
          sig_mean_mll, sig_sigmal_mll, sig_sigmar_mll,
          sig_alphal_mll, sig_alphar_mll, sig_nl_mll, sig_nr_mll]:
    print(f"  {v.GetName():30s} = {v.getVal():.6f}")

# ─────────────────────────────────────────────────────────────────────────────
# Component PDFs  (mirroring fit2D.py)
# ─────────────────────────────────────────────────────────────────────────────
bkgnonpeak_met = ROOT.RooGenericPdf(
    "bkgnonpeak_met", "bkgnonpeak_met",
    "1/b_nonpeak_met * exp(-(@0 - mu_nonpeak_met)/b_nonpeak_met"
    " - exp(-(@0 - mu_nonpeak_met)/b_nonpeak_met))",
    ROOT.RooArgList(met, mu_nonpeak_met, b_nonpeak_met))

bkgnonpeak_mll = ROOT.RooExponential("bkgnonpeak_mll", "bkgnonpeak_mll",
                                      mll, a_nonpeak_mll)

bkgpeaking_met = ROOT.RooGenericPdf(
    "bkgpeaking_met", "bkgpeaking_met",
    "1/b_peaking_met * exp(-(@0 - mu_peaking_met)/b_peaking_met"
    " - exp(-(@0 - mu_peaking_met)/b_peaking_met))",
    ROOT.RooArgList(met, mu_peaking_met, b_peaking_met))

bkgpeaking_mll = ROOT.RooCrystalBall("bkgpeaking_mll", "bkgpeaking_mll",
                                      mll, zp_mean_mll, zp_sigmal_mll, zp_sigmar_mll,
                                      zp_alphal_mll, zp_nl_mll, zp_alphar_mll, zp_nr_mll)

sig_dcb_mll = ROOT.RooCrystalBall("sig_dcb_mll", "sig_dcb_mll",
                                   mll, sig_mean_mll, sig_sigmal_mll, sig_sigmar_mll,
                                   sig_alphal_mll, sig_nl_mll, sig_alphar_mll, sig_nr_mll)

pdf_of_spline = workspace.pdf("pdf_of_spline_650_1")

# ─────────────────────────────────────────────────────────────────────────────
# 2D PDFs
# ─────────────────────────────────────────────────────────────────────────────
bkgnonpeak_2dpdf = ROOT.RooProdPdf("bkgnonpeak_2dpdf", "bkgnonpeak_2dpdf",
                                    ROOT.RooArgList(bkgnonpeak_mll, bkgnonpeak_met))
bkgpeaking_2dpdf = ROOT.RooProdPdf("bkgpeaking_2dpdf", "bkgpeaking_2dpdf",
                                    ROOT.RooArgList(bkgpeaking_mll, bkgpeaking_met))
sig_2dpdf        = ROOT.RooProdPdf("sig_2dpdf", "sig_2dpdf",
                                    ROOT.RooArgList(sig_dcb_mll, pdf_of_spline))

# ─────────────────────────────────────────────────────────────────────────────
# Floating parameters and extended S+B model
# r is fixed at TRUE_R; only n_sig and n_bkg float in the fit.
# n_sig and n_bkg initial values are the true yields used for generation.
# ─────────────────────────────────────────────────────────────────────────────
r     = ROOT.RooRealVar("r",     "Peaking bkg fraction",  TRUE_R,    0.0,  1.0)
n_sig = ROOT.RooRealVar("n_sig", "Signal yield",          TRUE_N_SIG, -5.0, 10.0)
n_bkg = ROOT.RooRealVar("n_bkg", "Background yield",      TRUE_N_BKG,  0.0, 50.0)

r.setConstant(True)

bkgtot_2dpdf = ROOT.RooAddPdf("bkgtot_2dpdf", "bkgtot_2dpdf",
                                ROOT.RooArgList(bkgpeaking_2dpdf, bkgnonpeak_2dpdf), r)

model = ROOT.RooAddPdf("model", "model",
                        ROOT.RooArgList(sig_2dpdf, bkgtot_2dpdf),
                        ROOT.RooArgList(n_sig, n_bkg))

# ─────────────────────────────────────────────────────────────────────────────
# RooMCStudy: generate and fit
# generateAndFit(nSamples, nEvtPerSample=0, keepGenData=True)
#   nEvtPerSample=0 → extended mode: Poisson(n_sig + n_bkg) events per toy
# ─────────────────────────────────────────────────────────────────────────────
ROOT.RooRandom.randomGenerator().SetSeed(args.seed)

mcs = ROOT.RooMCStudy(
    model, obs_set,
    RF.Extended(True),
    RF.Silence(),
    RF.FitOptions(RF.Extended(True), RF.PrintLevel(-1), RF.Strategy(1)),
)

print(f"Running generateAndFit({args.n_toys}) with keepGenData=True ...")
mcs.generateAndFit(args.n_toys, 0, True)
print("Done.")

# ─────────────────────────────────────────────────────────────────────────────
# Per-parameter binnings for Canvas 1
#   Bin widths: n_sig param=0.5, n_bkg param=1, n_sig pull=1, n_bkg pull=0.2
# ─────────────────────────────────────────────────────────────────────────────
_nsig_vals, _nbkg_vals = [], []
for _i in range(args.n_toys):
    _fr = mcs.fitResult(_i)
    if _fr and _fr.status() == 0:
        _p = _fr.floatParsFinal()
        _nsig_vals.append(_p.find("n_sig").getVal())
        _nbkg_vals.append(_p.find("n_bkg").getVal())

def _auto_binning(vals, bw):
    lo = bw * math.floor((min(vals) - bw) / bw)
    hi = bw * math.ceil((max(vals) + bw) / bw)
    return ROOT.RooBinning(int(round((hi - lo) / bw)), lo, hi)

_fallback_nsig = ROOT.RooBinning(20, n_sig.getMin(), n_sig.getMax())
_fallback_nbkg = ROOT.RooBinning(20, 0.0, 40.0)

param_binning = {
    "n_sig": _auto_binning(_nsig_vals, 0.5) if _nsig_vals else _fallback_nsig,
    "n_bkg": _auto_binning(_nbkg_vals, 1.0) if _nbkg_vals else _fallback_nbkg,
}
pull_binning = {
    "n_sig": ROOT.RooBinning(int(round(8 / 1.0)), -4.0, 4.0),
    "n_bkg": ROOT.RooBinning(int(round(8 / 0.2)), -4.0, 4.0),
}

# ─────────────────────────────────────────────────────────────────────────────
# Plotting helpers
# ─────────────────────────────────────────────────────────────────────────────
keep_alive = []

TITLE_SIZE = 0.075
LABEL_SIZE = 0.065
Y_OFFSET   = 1.30

def style_frame(frame):
    frame.GetXaxis().SetTitleSize(TITLE_SIZE)
    frame.GetXaxis().SetLabelSize(LABEL_SIZE)
    frame.GetYaxis().SetTitleSize(TITLE_SIZE)
    frame.GetYaxis().SetLabelSize(LABEL_SIZE)
    frame.GetYaxis().SetTitleOffset(Y_OFFSET)

# ─────────────────────────────────────────────────────────────────────────────
# Canvas 1 — fitted parameter distributions + pulls  (2×2)
#   Top row:    n_sig distribution,  n_bkg distribution
#   Bottom row: n_sig pull,          n_bkg pull
# ─────────────────────────────────────────────────────────────────────────────
c_mcs = ROOT.TCanvas(f"toy_mcs_{args.mode}", f"Toy study ({args.mode})", 1200, 600)
c_mcs.Divide(2, 2, 0.002, 0.002)

true_vals = {n_sig.GetName(): TRUE_N_SIG, n_bkg.GetName(): TRUE_N_BKG}

for col, var in enumerate([n_sig, n_bkg], start=1):
    # ── fitted parameter distribution (row 1) ──────────────────────────────
    c_mcs.cd(col)
    ROOT.gPad.SetLeftMargin(0.18)
    ROOT.gPad.SetBottomMargin(0.18)
    frame_param = mcs.plotParam(var, RF.Binning(param_binning[var.GetName()]), RF.FitGauss(True))
    style_frame(frame_param)
    frame_param.Draw()
    tval = true_vals[var.GetName()]
    true_line = ROOT.TLine(tval, 0, tval, frame_param.GetMaximum() * 0.9)
    true_line.SetLineColor(ROOT.kRed + 1)
    true_line.SetLineStyle(ROOT.kDashed)
    true_line.SetLineWidth(2)
    true_line.Draw("SAME")
    keep_alive.extend([frame_param, true_line])

    # ── pull distribution (row 2) ───────────────────────────────────────────
    c_mcs.cd(col + 2)
    ROOT.gPad.SetLeftMargin(0.18)
    ROOT.gPad.SetBottomMargin(0.18)
    frame_pull = mcs.plotPull(var, RF.Binning(pull_binning[var.GetName()]), RF.FitGauss(True))
    style_frame(frame_pull)
    frame_pull.Draw()
    zero_line = ROOT.TLine(0, 0, 0, frame_pull.GetMaximum() * 0.9)
    zero_line.SetLineColor(ROOT.kBlack)
    zero_line.SetLineStyle(ROOT.kDashed)
    zero_line.SetLineWidth(2)
    zero_line.Draw("SAME")
    keep_alive.extend([frame_pull, zero_line])

c_mcs.cd(0)
info_mcs = ROOT.TLatex()
info_mcs.SetNDC()
info_mcs.SetTextSize(0.030)
info_mcs.SetTextFont(42)
info_mcs.SetTextAlign(22)
info_mcs.DrawLatex(
    0.50, 0.988,
    f"{args.mode.upper()} | {args.n_toys} toys | "
    f"true: n_{{sig}}={TRUE_N_SIG}, n_{{bkg}}={TRUE_N_BKG}, r={TRUE_R} (fixed)")
keep_alive.append(info_mcs)

outname = f"toy_study_{args.mode}_{args.n_toys}toys"
c_mcs.SaveAs(f"{outname}.png")
c_mcs.SaveAs(f"{outname}.pdf")

# ─────────────────────────────────────────────────────────────────────────────
# Canvas 3 — box-and-whisker plot of per-bin event counts across all toy datasets
#   x-axis = observable bin center  (m(ll) or MET)
#   y-axis = event count in that bin
#   Each box (one per observable bin) shows median, IQR, and whiskers over toys.
#   Red markers show the expected per-bin count at the true yields.
# ─────────────────────────────────────────────────────────────────────────────
n_generated = args.n_toys
print(f"Building violin plot from {n_generated} generated datasets ...")

N_BINS_TOY_MLL = mll.getBins()
N_BINS_TOY_MET = met.getBins()

toy_hmll, toy_hmet = [], []
for i in range(n_generated):
    ds_i = mcs.genData(i)
    h_mll_i = ROOT.TH1F(f"hmll_toy{i}", "", N_BINS_TOY_MLL, mll.getMin(), mll.getMax())
    h_met_i = ROOT.TH1F(f"hmet_toy{i}", "", N_BINS_TOY_MET, met.getMin(), met.getMax())
    h_mll_i.Sumw2(False)
    h_met_i.Sumw2(False)
    for j in range(ds_i.numEntries()):
        row = ds_i.get(j)
        h_mll_i.Fill(row.find(mll.GetName()).getVal())
        h_met_i.Fill(row.find(met.GetName()).getVal())
    toy_hmll.append(h_mll_i)
    toy_hmet.append(h_met_i)

# y-axis range: 0 … max observed count + a little headroom
max_count_mll = int(max(h.GetMaximum() for h in toy_hmll))
max_count_met = int(max(h.GetMaximum() for h in toy_hmet))

# TH2: x = observable bin center, y-column = distribution of counts across toys.
# ROOT draws each x-bin column as one violin with Draw("VIOLIN").
th2_vln_mll = ROOT.TH2F("th2_vln_mll", "",
    N_BINS_TOY_MLL, mll.getMin(), mll.getMax(),
    max_count_mll + 2, -0.5, max_count_mll + 1.5)
th2_vln_met = ROOT.TH2F("th2_vln_met", "",
    N_BINS_TOY_MET, met.getMin(), met.getMax(),
    max_count_met + 2, -0.5, max_count_met + 1.5)

for h in toy_hmll:
    for b in range(1, N_BINS_TOY_MLL + 1):
        th2_vln_mll.Fill(h.GetBinCenter(b), h.GetBinContent(b))
for h in toy_hmet:
    for b in range(1, N_BINS_TOY_MET + 1):
        th2_vln_met.Fill(h.GetBinCenter(b), h.GetBinContent(b))

for th2 in [th2_vln_mll, th2_vln_met]:
    th2.SetFillColorAlpha(ROOT.kBlue + 1, 0.45)
    th2.SetLineColor(ROOT.kBlue + 1)
    th2.SetLineWidth(1)

# Expected per-bin count: large generated sample scaled to true yield
expected_total = TRUE_N_SIG + TRUE_N_BKG
ROOT.RooRandom.randomGenerator().SetSeed(args.seed + 99999)
N_EXPECT = 50000
ds_expect = model.generate(obs_set, N_EXPECT)
h_exp_mll = ROOT.TH1F("h_exp_mll", "", N_BINS_TOY_MLL, mll.getMin(), mll.getMax())
h_exp_met = ROOT.TH1F("h_exp_met", "", N_BINS_TOY_MET, met.getMin(), met.getMax())
h_exp_mll.Sumw2(False)
h_exp_met.Sumw2(False)
for j in range(ds_expect.numEntries()):
    row = ds_expect.get(j)
    h_exp_mll.Fill(row.find(mll.GetName()).getVal())
    h_exp_met.Fill(row.find(met.GetName()).getVal())
h_exp_mll.Scale(expected_total / N_EXPECT)
h_exp_met.Scale(expected_total / N_EXPECT)
for h in [h_exp_mll, h_exp_met]:
    h.SetLineColor(ROOT.kRed + 1)
    h.SetLineWidth(2)
    h.SetMarkerColor(ROOT.kRed + 1)
    h.SetMarkerStyle(ROOT.kFullCircle)
    h.SetMarkerSize(1.0)
    h.SetFillStyle(0)

c_toys = ROOT.TCanvas(f"toys_overlay_{args.mode}",
                      f"Toy overlay ({args.mode})", 1000, 500)
c_toys.Divide(2, 1, 0.002, 0.002)

for pad_idx, (th2, h_exp, xlabel) in enumerate(
        [(th2_vln_mll, h_exp_mll, "m(ll) [GeV]"),
         (th2_vln_met, h_exp_met, "MET [GeV]")], start=1):
    c_toys.cd(pad_idx)
    ROOT.gPad.SetLeftMargin(0.18)
    ROOT.gPad.SetBottomMargin(0.18)
    th2.GetXaxis().SetTitle(xlabel)
    th2.GetYaxis().SetTitle("Events / bin")
    th2.GetXaxis().SetTitleSize(TITLE_SIZE)
    th2.GetXaxis().SetLabelSize(LABEL_SIZE)
    th2.GetYaxis().SetTitleSize(TITLE_SIZE)
    th2.GetYaxis().SetLabelSize(LABEL_SIZE)
    th2.GetYaxis().SetTitleOffset(Y_OFFSET)
    th2.Draw("CANDLE")
    h_exp.Draw("P SAME")

c_toys.cd(0)
info_toys = ROOT.TLatex()
info_toys.SetNDC()
info_toys.SetTextSize(0.030)
info_toys.SetTextFont(42)
info_toys.SetTextAlign(22)
info_toys.DrawLatex(
    0.50, 0.988,
    f"{args.mode.upper()} | {n_generated} toys (box-and-whisker) | "
    f"expected {expected_total:.1f} events (red)")
keep_alive.extend([c_toys, th2_vln_mll, th2_vln_met, h_exp_mll, h_exp_met,
                   ds_expect, info_toys] + toy_hmll + toy_hmet)

toys_outname = f"toy_study_overlay_{args.mode}_{args.n_toys}toys"
c_toys.SaveAs(f"{toys_outname}.png")
c_toys.SaveAs(f"{toys_outname}.pdf")

eosdir = "/eos/user/s/skkwan/www/higgsino/studies/mll-MET-fit-2D/toys"
os.system(f"mv {outname}.* {eosdir}")
os.system(f"mv {toys_outname}.* {eosdir}")
print(f"Moved outputs to {eosdir}")
