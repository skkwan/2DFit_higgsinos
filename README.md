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
Check the instructions in `hbb_zll`.