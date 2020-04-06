# ODUSSEAS (Observing Dwarfs Using Stellar Spectroscopic Energy-Absorption Shapes):
# A Machine-Learning Tool for the derivation of Teff and [Fe/H] of M Dwarf stars

IF YOU USE THIS TOOL IN YOUR RESEARCH, PLEASE CITE THE CORRESPONDING PAPER:

https://www.aanda.org/articles/aa/abs/2020/04/aa37194-19/aa37194-19.html

"ODUSSEAS.py" is the code we run.

We can set the regression type using the setting "regression".
This can be: 'ridge' (recommended), 'ridgecv', 'linear', 'multitasklasso', 'multitaskelasticnet'
We can also choose to do r.v. correction to our spectra if they are shifted, by setting to 'yes' the "do_rv_cor" option.

Input: inside a folder with the path "spectra/newstars/", there should be the fits files of the 1D spectra of the unknown stars. Their filepaths should be written in a text file called "1Dfilelist.dat", and next to them the resolution of each spectrum. See example below:

spectra/newstars/starA.fits 115000

spectra/newstars/starB.fits 94600

spectra/newstars/starC.fits 75000

Output: A text file named "Parameter_Results.dat" is created. It contains the averages of [Fe/H] and Teff after 100 M.L.runs for each star along with their dispersion, the scores and the mean absolute errors of the models that predicted them.

Demo set: 1D spectra of stars from 5 different spectrographs with different resolutions and respective HARPS datasets for them are provided to use our tool. For comparison, the reference values of the respective HARPS spectra are the following: Gl846 = -0.08 & 3682; Gl514 = -0.13 & 3574 ; Gl908 = -0.38 & 3587 ; Gl674 = -0.18 & 3284 and for the HARPS star outside the reference HARPS dataset Gl643 = -0.26 & 3102 by Neves et al (2014).

We already provide precomputed pseudo EWs for a range of spectral resolutions used in popular spectrographs. For completeness, the repository also includes the code "HARPS_dataset.py", which can create a library of M dwarfs from our HARPS sample for any resolution we want to work at (the associated fits files are not uploaded). If you wish to create additional libraries, please contact us.

Note: The structure of the repository and associated files should be kept as found, in order for the code to run properly.

