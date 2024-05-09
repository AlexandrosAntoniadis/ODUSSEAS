import MachineLearning
import New_data

reference = "interferometry"  # choose the reference scale: 'photometry' for 65 stars with Teff from Casagrande08 and [Fe/H] from Neves12, or 'interferometry' for 47 stars with Teff from Khata21 and Rabus19, and [Fe/H] from Neves12
regression = "ridge"  # choose the ML model: linear, ridge, ridgecv, multitasklasso, multitaskelasticnet  (Recommended: ridge)
do_rv_cor = "yes"  # choose 'yes' if you want to do rv correction to the spectra or 'no' if they are already corrected

New_data.EWmeasurements(do_rv_cor)
MachineLearning.ML(regression, reference)
