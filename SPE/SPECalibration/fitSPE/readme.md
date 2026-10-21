These three files are used to fit the single-photoelectron (SPE) amplitude distributions from SBND September 2026 run 20530, using a single Gaussian for each PMT channel. They work with the ROOT histogram file produced by the PMT gain calibration module.

|File|Purpose|
|---|---|
|`fitspe_sep.C`|Automatically fits SPE amplitude histograms for all channels found in the input ROOT file and produces diagnostic plots and a results table.|
|`fitspe_manual_sep.C`|Manually refits selected channels after reviewing the automatic results, using user-defined fitting windows.|
|`sep_voltage_settings.h`|Contains the September nominal operating voltage for each PMT, its target SPE amplitude (12.75 or 14.00 ADC), and off / hardware-voltage-limit flags. Both fitting macros include this header.|
Before running, **change the directory in both** `fitspe_sep.C` **and** `fitspe_manual_sep.C` to your own output directory:
```
TString base = gSystem->ExpandPathName("~/work/output/Sep/");
```
Place `pmt_gain20530.root` in that directory.

## 1. Automatic Gaussian fitting
```
root -l -b -q fitspe_sep.C
```

Review `fitspe_sep.pdf` to identify channels whose automatic fitting window or Gaussian fit needs adjustment. A flag is a reason to inspect the fit, not necessarily proof that it is unusable.

## 2. Manual refitting (if needed)
Open `fitspe_manual_sep.C` and add or edit entries in the `MANUAL` map for channels that require refitting:
```
static std::map<int,std::array<double,3>> MANUAL = {
    {140, {10.0, 18.0, 0.0}},
};
```
Each entry has the form:
```
{channel, {fit_lower_bound_ADC, fit_upper_bound_ADC, initial_mean_ADC}}
```
Here, an `initial_mean` of `0.0` lets the macro choose the strongest histogram bin within the specified fitting interval as its starting mean.

Run:
```
root -l -b -q fitspe_manual_sep.C
```
This creates `fitspe_sep_manual.pdf` with plots of the selected manual fits and updates the corresponding channel rows in `fitspe_sep.txt` in place. Other channel rows are retained. The manual script requires the automatic fitting results file to exist first.

## September voltage settings
`sep_voltage_settings.h` provides the per-channel settings used for September run 20530:
- `voltage`: nominal operating HV in volts (integer value used for the run).
- `target`: the intended SPE amplitude, **12.75 or 14.00 ADC** for operating PMTs.
- `off`: whether the PMT was switched off (listed with zero voltage).
- `atLimit`: whether the HV is at the PMT hardware limit (1780 or 1790 V in this table). Limit-only entries are assigned a 12.75 ADC target in the header.