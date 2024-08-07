# ODUSSEAS (Observing Dwarfs Using Stellar Spectroscopic Energy-Absorption Shapes):
## A Machine-Learning Tool for the derivation of Teff and [Fe/H] of M Dwarf stars

IF YOU USE THIS TOOL IN YOUR RESEARCH, PLEASE CITE THE CORRESPONDING PAPER:

https://doi.org/10.1051/0004-6361/201937194


## Usage
```bash
$ ODUSSEAS --help

 Usage: ODUSSEAS [OPTIONS] INPUT_SPECTRA

 Run ODUSSEAS with the arguments as listed below

╭─ Arguments ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ *    input_spectra      TEXT  [default: None] [required]                                                                                         │
╰──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
╭─ Options ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ --reference                                            [photometry|interferometry]                  choose the reference scale: 'photometry' for │
│                                                                                                     65 stars with Teff from Casagrande08 and     │
│                                                                                                     [Fe/H] from Neves12, or 'interferometry' for │
│                                                                                                     47 stars with Teff from Khata21 and Rabus19, │
│                                                                                                     and [Fe/H] from Neves12                      │
│                                                                                                     [default: interferometry]                    │
│ --regression                                           [linear|ridge|ridgecv|multitasklasso|multit  choose the ML model. Recommended: ridge      │
│                                                        askelasticnet ]                              [default: ridge]                             │
│ --verbose                 --no-verbose                                                              [default: no-verbose]                        │
│ --skip-ew-measurements    --no-skip-ew-measurements                                                 If this step is already done, then it can be │
│                                                                                                     skipped in further analysis, as it is a bit  │
│                                                                                                     slow                                         │
│                                                                                                     [default: no-skip-ew-measurements]           │
│ --help                                                                                              Show this message and exit.                  │
╰──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
```

### Example
```bash
$ ODUSSEAS 1Dfilelist.dat
```
An example of the `1Dfilelist.dat` file can be seen in this repository.

This should run the example provided with the code. It should output the
results fom 5 spectra in the folders `results/` (the pseudo EWs measured for
each spectrum), and `Parameter_Results.dat` (the calculated atmospheric
parameters for each star).

## Installation
It is recommended to install this package in a virtual or conda environment
following with the command
```bash
$ pip install git+https://github.com/AlexandrosAntoniadis/ODUSSEAS
```
A recent version of python should work, but do let us now if you have any
issues installing and running the code.

You can also clone this repository and install it locally.
```bash
$ git clone https://github.com/AlexandrosAntoniadis/ODUSSEAS
$ cd ODUSSEAS
$ pip install -e .
```


## Documentation
We select the methods by which the reference parameters have been derived,
using the setting `reference`. This can be: `photometry` which uses as
reference dataset 65 stars with photometric scales of Teff by Casagrande et al.
(2008) and [Fe/H] by Neves et al. (2012), or `interferometry` (regarded as the
new version of ODUSSEAS) which uses as reference dataset 47 stars with
interferometry-based Teff by Khata et al. (2021) and Rabus et al. (2019) and
[Fe/H] derived with the method by Neves et al. (2012) using the updated
parallaxes from Gaia DR3. We can set the regression type using the setting
`regression`. This can be: `ridge` (recommended), `ridgecv`, `linear`,
`multitasklasso`, `multitaskelasticnet`

Input: inside a folder with the path "spectra/newstars/", there should be the
fits files of the 1D spectra of the unknown stars. Their filepaths should be
written in a text in same format as `1Dfilelist.dat`, and next to them the
resolution of each spectrum. See example below:

```
spectra/newstars/starA.fits 115000
spectra/newstars/starB.fits 94600
spectra/newstars/starC.fits 75000
```

Output: A text file named `Parameter_Results.dat` is created. It contains the
average values of [Fe/H] and Teff after 100 M.L. runs for each star, along with
their dispersion, the mean absolute errors of the models that predicted them,
the wide error budget (after taking into consideration the intrinsic
uncertainties of the reference parameters into the machine learning process),
and the machine-learning scores.

Demo set: 1D spectra of stars from 5 different spectrographs with different
resolutions and respective HARPS datasets for them are provided to use our
tool.
For comparison, the reference values of the respective HARPS spectra are the following:
Using the scales of Casagrande08 and Neves12: Gl846 = -0.08 & 3682 ; Gl514 =
-0.13 & 3574 ; Gl908 = -0.38 & 3587 ; Gl674 = -0.18 & 3284 and for the HARPS
star outside the reference HARPS dataset Gl643 = -0.26 & 3102 by Neves et al
(2014).
Using the scales of Khata21 & Rabus19 and updated Neves12: Gl846 = -0.07 & 3810
; Gl514 = -0.15 & 3671 ; Gl908 = -0.40 & 3475 ; Gl674 = -0.19 & 3409 ; Gl643 =
-0.32 & 3243.

We already provide precomputed pseudo EWs for a range of spectral resolutions
used in popular spectrographs. For completeness, the repository also includes
the code `src/odusseas/examples/HARPS_dataset.py`, which can create a library
of M dwarfs from our HARPS sample for any resolution we want to work at (the
associated fits files are not uploaded). If you wish to create additional
libraries or for any other question, please contact us at:
alexandros.antoniadis@astro.up.pt or open an issue on Github.
