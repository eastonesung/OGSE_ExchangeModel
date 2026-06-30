# OGSE Exchange-Sensitized Kurtosis Imaging

MATLAB code accompanying:

> Sung D, et al. **Microstructural and exchange imaging with oscillating gradient
> spin-echo (OGSE) diffusion MRI.** *Magnetic Resonance in Medicine*, 2025.
> [https://doi.org/10.1002/mrm.70300](https://doi.org/10.1002/mrm.70300)

This repository contains the model, fitting, and figure-generation code used in
the paper. It is shared for transparency and reproducibility; it is **not**
a turnkey package, since it depends on lab-internal data paths and utility
code (see [Dependencies](#dependencies) below) that you will need to
configure for your own environment.

## Repository structure

```
.
├── preprocessing/
    ├── pipeline_ogse.sh                       Entire pipeline including preprocessing and dti/dki fitting
    ├── s01_rung_preprocessing_ogse_C2.sh      Preprocessing based on DESIGNER pipeline
    └── s02_dti_dki_analysis.sh                DTI/DKI fitting using MRtrix3 as well as FSL

├── models/
│   └── gpu_KurtModel.m                       Forward model class (askAdam/MCMC compatible):
│                                              frequency-dependent kurtosis as a function of
│                                              exchange time, volume fraction, tortuosity, etc.
├── analysis/
│   └── roi_kurtosis_fitting_askAdam_MCMC.m   ROI-level extraction of mean kurtosis from
│                                              FreeSurfer/SynthSeg parcellations across the
│                                              OGSE frequency protocol, and model fitting.
└── figures/
    ├── make_manuscript_figures.m             Main + supplementary manuscript figures
    │                                          (waveforms, power spectra, Kvar/h[q]/MK vs.
    │                                          frequency, mixing-time optimization, group MD/MK).
    └── make_waterphantom_figures.m           Water-phantom validation figures (MD/MK vs.
                                               frequency in an isotropic, no-exchange medium).
```

## Dependencies

These scripts are **not** standalone — install/clone the following first and
point to them in the "USER CONFIGURATION" block at the top of each script:

1. **[gacelle](https://github.com/kschan0214/gacelle)** (GPU-AcCELerated toolbox
   for high-throughput multi-dimensional quantitative parameter Estimation).
   Provides the `askadam` (gradient-descent) and `mcmc` solvers and the
   `utils` helper class used throughout this code.
2. **`util_ogse`** — an author-maintained OGSE waveform utility class
   (cosine-trapezoidal waveform generation, q(t)/q(w) calculation, `KMhq`,
   `KvarCosTrap`, `KMhqOG`, `plot_waveform`, etc.). This class is referenced
   throughout the scripts (`ut = util_ogse()`) but **is not yet included in
   this repository** — add it under e.g. `util_ogse/` and point
   `util_ogse_dir` at it, or fold its contents into the `gacelle` install.
3. MATLAB toolboxes: **Image Processing** (`niftiread`/`niftiinfo`),
   **Parallel Computing** (`gpuArray`, `parfor`), and **Deep Learning**
   (`dlarray`) for GPU-accelerated model fitting.

> If you (the repo owner) have `util_ogse` and a minimal set of de-identified
> example data, consider adding both — without them, users can read the code
> but cannot run it end-to-end.

## Data

The scripts expect a BIDS-organized derivatives tree (`bids/derivatives/...`)
with preprocessed DWI, FreeSurfer/SynthSeg segmentations, and MRtrix3 tensor
outputs for each subject. Subject IDs (e.g. `sub-ogse006`) are anonymized
study identifiers. Raw/derivative imaging data are not included in this
repository; see the data availability statement in the paper for access
information.

## Usage

1. Edit the `USER CONFIGURATION` block at the top of each script
   (`gacelle_dir`, `util_ogse_dir`, `project_dir`).
2. Run `analysis/roi_kurtosis_fitting_askAdam_MCMC.m` to extract ROI kurtosis
   and fit the exchange model.
3. Run the scripts in `figures/` (cell-by-cell — each `%%` section
   corresponds to one figure/panel) to reproduce the manuscript figures.

## Citation

If you use this code, please cite the paper above, and please also cite
**gacelle** per the instructions in its
[repository](https://github.com/kschan0214/gacelle#terms-of-use), since the
fitting routines here depend on it.

## License

Released under the MIT License (see `LICENSE`). Note that the `gacelle`
dependency is licensed separately under GPL-3.0 — see that repository for
its terms.

## Contact

Dongsuk Sung — dsung2@mgh.harvard.edu
Athinoula A. Martinos Center for Biomedical Imaging, MGH/Harvard Medical School
