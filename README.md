# ODUSSEAS (Observing Dwarfs Using Stellar Spectroscopic Energy-Absorption Shapes):
# A Machine-Learning Tool for the derivation of Teff and [Fe/H] of M Dwarf stars

"ODUSSEAS.py" is the code we run.
In its top, there is a setting called "option".

It can be either:
"measure" - the code both measures the pseudo EWs and derives the stellar parameters of unknown stars by machine learning ;

or:
"justML" - the code derives the stellar parameters by machine learning (requires previously measured pseudo EWs). Use for quick runs, useful to produce new models applied to each star, leads to new calculations of parameters. To be used when the first run does not return very high scores for a star measurement, therefore used for producing more precise models and parameter values.

We can also set the regression type using the setting "regression".
This can be: "ridge" (recommended), "ridgecv", "linear", "multitasklasso", "multitaskelasticnet"

Input: inside a folder with the path "spectra/newstars/", there should be the fits files of the 1D spectra of the unknown stars (corrected for R.V.). Their filepaths should be written in a text file called "1Dfilelist.dat", and next to them the resolution of each spectrum. See example below:

spectra/newstars/starA.fits 115000

spectra/newstars/starB.fits 94600

spectra/newstars/starC.fits 75000

Output: A text file named "Parameter_Results.dat" is created. It contains the [Fe/H] and Teff of the unknown stars, along with the scores and mean absolute errors of the models that predicted each one of them.

Demo set: 1D spectra of stars from 5 different spectrographs with different resolutions and respective HARPS datasets for them is provided to use our tool. For comparison, the reference values (obtained from their HARPS spectra) are the following: Gl402 = 0.06 & 2984 ; GJ163 = 0.06 & 3274 ; GJ1061 = -0.02 & 2938 ; Gl514 = -0.16 & 3539 ; Gl408 = -0.47 & 3530 and for the HARPS star outside the HARPS dataset GJ3470 = 0.08 & 3470 by Santos et al (2013).

We already provide precomputed pseudo EWs for a range of spectral resolutions used in popular spectrographs. For completeness, the repository also includes the code "HARPS_dataset.py", which can create a library of M dwarfs from our HARPS sample for any resolution we want to work at (the associated fits files are not uploaded). If you wish to create additional libraries, please contact us.

Note: The structure of the repository and associated files should be kept as found, in order for the code to run properly.
