These ROOT macros fit the pedestal and SPE peaks for the August 2026 voltage scan.

## Fitting models

| Macro | Combined model | Pedestal fit region |
| --- | --- | --- |
| `combfit_gaussian_aug.C` | Gaussian pedestal + Gaussian SPE | Right side of pedestal, starting at its measured peak |
| `combfit_expo_aug.C` | Exponential pedestal + Gaussian SPE | Right side of pedestal, starting at its measured peak |
| `combfit_gamma_aug.C` | Shifted-gamma pedestal + Gaussian SPE | Includes the rising side and the peak of the pedestal |

## Inputs and paths

All three macros read `PMTGain/amp_opchannel_*` histograms from:

- `pmt_gain20468.root` — new nominal HV

- `pmt_gain20469.root` — nominal +10 V

- `pmt_gain20470.root` — nominal +20 V

- `pmt_gain20471.root` — nominal +30 V

The default input directory is `~/work/output/Aug_scan/`. **Edit `base` and `outDir` near the start of each macro** to match your directories. Outputs are written to `~/work/output/combfit_aug/` by default.

If `fitspe_aug.txt` is available in the input directory, the macros use its single-Gaussian results as SPE initial values and for comparison plots. Otherwise, they estimate SPE seeds from the histograms; the saved single-fit comparison curves may be unavailable.

## Run

From the directory containing the `.C` files, execute each macro separately:

```bash

root -l -b -q 'combfit_gaussian_aug.C'

root -l -b -q 'combfit_expo_aug.C'

root -l -b -q 'combfit_gamma_aug.C'

```

## Outputs

Each macro writes four files, replacing `<model>` with `gaussian`, `expo`, or `gamma`:

- `combfit_<model>_aug.txt` — fitted pedestal/SPE parameters, SPE mean and error, fit status, quality flags, and peak-to-valley metrics.

- `combfit_<model>_aug.pdf` — fitted distributions for the four voltages, grouped by PMT channel.

- `combfit_<model>_aug_residual.pdf` — fits and standardized residuals.

- `combfit_<model>_aug_vs_single.pdf` — comparison with the saved single-Gaussian SPE fits when available.

The `CLEAR` / `BURIED` labels are based on an observed peak-to-valley ratio (default threshold `PV_RATIO_MIN = 1.10`); they are diagnostic labels, not a manually validated exclusion list.