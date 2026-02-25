# 2D fit repository

Original code from [upstream by Meraj](https://github.com/mhussainphys/2DFit_higgsinos).

## Status
- Only 2018 `mm` channel possible at the moment because `weight_nominal_mm` and `weight_nominal_ee` are named differently and cannot be used at the same time in the current script
- Fit parameter result values are hard-coded into plotting 

## To run the original m(gamma gamma) + MET fit
1. In the top-level directory:
```bash
python3 py_2d_SonlyFit_BonlyFit.py
```

## To work on the m(ll) + MET fit
1. `cd hbb_zll/`
2. Copy the input signal file to the current directory (e.g. `snapshot_TChiZH_650_1_cat_0_batch_0_channel_mm_SR_mll_MET_fit_scheme.root`): currently only `mm` channel is supported
3. Hadd the background MC files:
```bash
# Point to the correct directory in here
python3 reformat.py
```
4. Run the fit:
```bash
python3 py_2d_fit.py
```
5. Currently the fit parameter results are hard-coded into the plotting:
```bash
python3 plot1DProjection.py
```