# ODUSSEAS (Observing Dwarfs Using Stellar Spectroscopic Energy-Absorption Shapes):
# A Machine-Learning Tool for the derivation of Teff and [Fe/H] of M Dwarf stars

"ODUSSEAS.py" makes use of two algorithms: the "Newdataset.py", for the measurement the pseudo EWs of unknown spectra, and the "MachineLearning.py" for the derivation of their [Fe/H] and effective temperature. 

Inside the folder called "spectra", the user creates a second folder called "newstars", where he keeps the fits files of the 1D spectra of the unknown stars (having already corrected for radial velocity). These filepaths shall be written in the file called "1Dfilelist.dat" and next to them the resolution of each star.

In the top of the code "ODUSSEAS.py", there are three settings: 

1st) The main setting  is called "option". It shall be setted as "measure" if the user want to measure the new spectra for the first time, so the code firstly calculates their pseudo EWs and then derives the parameters by Machine Learning. 
It can also be setted as "justML" in the case that the user wants to quickly try another run of Machine Learning. This action will produce new models applied for each star, which lead to new calculations of parameters. This can be done in the case that the first run does not return very high scores for a star measurement, so the user can try instantly another run that will produce a more precise model with more precise parameter values.

2nd) Furthermore, the user can set the resolution limit, considered as the convolution limit above which we convolve the new spectra to their own resolution. This is suggested in the case of spectra having very similar resolutions to the 115000 of the HARPS data set, since the act of convolution itself changes the shape of the absorption lines. We suggest this resolution limit to be 100000. 

3rd) Finally, the user can change the regression type to be used by the Machine Learning process. The "ridge" is recommended, but also "ridgecv" and "linear" work at similar level of efficiency as well. 

The workflow of "ODUSSEAS.py" goes as following: 
The "New_dataset.py" reads the filepaths and resolutions in the "1Dfilelist.dat" and it proceeds to either convolve the unknown spectra to their own resolution based on the convolution limit or not. The final spectra -either convolved or not- are created in the same folder, while automatically the "final1Dfilelist.dat" is created and read by the code. The code reads the wavelength range of each star individually and proceeds to the measurement of the respective pseudo-EWs which are saved as individual csv files for each star. These csv files are used during the operation of the "MachineLearning.py" that returns the values of [Fe/H] and Teff, saved in a file called "Parameter_Results.dat". 

***For testing the above, we provide 1D spectra from 5 different spectrographs with different resolutions and the respective HARPS csv files for these resolutions, based on which the Machine Learning will predict the parameters. 
For comparison, the reference values according to their HARPS spectra, are the following ( [Fe/H] & Teff] ): 
Gl402 = 0.06 & 2984 ;
GJ163 = 0.06 & 3274 ;
GJ3470 = 0.08 & 3470 ;
GJ1061 = -0.02 & 2938 ;
Gl514 = -0.16 & 3539 ;
Gl408 = -0.47 & 3530 ;


The "HARPSdataset.py" is the code that creates the library of the M dwarfs from the HARPS sample, for any resolution we want to work at.
It makes use of 96 stars with high-quality spectra from HARPS, in order to build machine-learning models as accurate and precise as possible, upon their measurements. Each time this code runs, the outcome is a csv file that contains the names of the stars, the absorption lines measured according to a linelist of 4104 lines from 530 to 690 nm, their pseudo EW values according to the resolution we convolve the spectra and their values of [Fe/H] and effective temperature that we got from calibration as reference values. 

At the top of the code, the settings the user does are simply two: The choice either to convolve the spectra to a new resolution or not (in the case that the csv file of this resolution does not already exists) and the resolution for the function of convolution to create the respective csv file. This output csv file of the resolution we want, is used later as input to the Machine Learning algorithm when running the Mdwarfs.py code for training the machine and testing the generated model.

More in detail, the procedure inside the "HARPSdataset.py" code goes as following: In the main folder called "spectra", we keep the "HARPS" folder which containes all the HARPS fits files (are not uploaded here). When the code runs, those 96 files of spectra with SNR > 25, are read with a "goodHARPSfilelist.dat" in order to construct their spectra based on their values of CRVAL1, CDELT1,  NAXIS1. They are convolved to the new resolution and the new fits files are saved with the respective names in the same folder. Automatically, the new filepaths are created and read for getting the spectra after convolution. The function "pseudoEW" starts to calculate all the pseudo EWs based on the linelist called "lines.rdb" that is read. The results of pseudo EWs for each star are saved in individual files inside a directory that is created automatically. Then, those values form a csv table which has as first column the names of the stars and first row the numbers of the central absorption lines (calculated as the middle by the left and right boundaries provided by the "lines.rdb" linelist). The final two columns adjusted to this csv table are the [Fe/H] and effective temperature values of each star, read by the file "originalparameters.dat", as they have been calculated by calibration and they are used as the reference values during the machine learning process.
