# Steps for fitting

1. **Initial Z-peak shape fit** (peaking background in m(ll), m(ll) shape): Combine the n-tuples used to perform the Z-peak fit:
```bash
# from the main directory
cd zpeak_fit/
python3 reformat_zPeak.py # creates backgrounds_CRZ_Zpeak_2018.root
python3 initial_ZPeak_fit.py # creates initial_zPeak_fit_result.root 
python3 plotZPeakFit.py # validate the results in a plot
```
From counting yields (i.e. number of events in the n-tuple weighed by the cross-section), about 0.0876 of the total background comes from events peaking in m(ll). 
Output plots are in [https://skkwan.web.cern.ch/higgsino/studies/mll-MET-fit-2D/zpeak_fit/](here).

2. **Get remaining background MET shapes (peaking and non-peaking MET shape, non-peaking m(ll) shape)**: We obtain the remaining three background PDFs with:
```bash
# from the main directory
cd individual_pdf_fits/
python3 reformat.py # creates backgrounds_peaking.root and backgrounds_nonpeak.root
python3 obtain_background_components.py # creates fitresult_background_all_except_ZPeak.root
```
Output plots are in [https://skkwan.web.cern.ch/higgsino/studies/mll-MET-fit-2D/background_shapes/](here).

3. **Get signal PDFs**:
```bash
# from the main directory
cd individual_pdf_fits/
python3 obtain_signal_components.py # creates fitresult_signal.root
```


3. 2fit to background, 2D fit to signal
Perform the 2D fit to background (using the Z-peak as the real mll component of the peaking-in-mll-background), and the 2D fit to signal:

```bash
# back in the main directory
cd overall_fit/
python3 reformat.py # prepare the input files: creates backgrounds_for_2D_fit.root and snapshot_TChiZH_650_1_SR_mll_MET_fit_scheme.root
python3 py_2d_fit_withZPeakFixed.py # creates fitresult.root
python3 plot1DProjection.py # plots the fit result and fit inputs 
```

## Next steps
- Background-only fit with S+B shapes, with the ratio of peaking and non-peaking backgrounds floating 
- Background+signal fit with S+B shapes, with the ratio of peaking and non-peaking backgrounds floating