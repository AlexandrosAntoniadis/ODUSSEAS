#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 31 12:54:57 2018

@author: mac
"""

import numpy as np
from astropy.io import fits
import matplotlib.pylab as plt
from scipy.interpolate import interp1d
from PyAstronomy import pyasl
from astropy.io import fits
from build_pew_verygeneral import pseudo_EW

import os
import pandas as pd


####################### 
# set the resolution that you want to convolve the HARPS ML sample.
# set the resonumber (the same as above in k units) to be added in the names of the files.
# define the instrument according to which, the same linelist range will be chosen.
#define the rangename as you want so to be added in the name of the final file

resolution = 110000
resonumber = "110"
linerange= "UVES"
#rangename= "smallerLR"

#######################

def convolve_data(fname, newresolution, resonumber):


    flux = fits.getdata(fname)
    hdr = fits.getheader(fname)
    w0, dw, N = hdr['CRVAL1'], hdr['CDELT1'], hdr['NAXIS1']
    wavelength = w0 + dw * np.arange(N)
    
    
    newflux = pyasl.instrBroadGaussFast(wavelength, flux, newresolution, edgeHandling="firstlast", fullout=False, maxsig=None)
    
    #prihdr = fits.Header()
    #prihdr.set("CDELT1",w0)
    #prihdr.set("CRVAL1",dw)
    #prihdr.set("NAXIS1",N)

    fits.writeto(fname.replace('.fits', '')+'res'+resonumber+'.fits', newflux, hdr, overwrite=True)

def Load_Data_Pairs(data1):
    
    # Get the data and the names from our results
    our_datanames = []
    our_data = []
    listdir = os.listdir(data1)
    listdir = np.sort(listdir)
    for file in listdir:
        # get all .dat files
        if file.endswith(".dat"):
            temp_name = os.path.join(data1, file)
            # remove unwanded elements and append the names
            our_datanames.append(temp_name.split('/')[1].strip( '.dat' ))
            # append data
            our_data.append(np.loadtxt(temp_name))
    
        
    return our_data, our_datanames   



filepaths = np.loadtxt('HARPSfilelist.dat', dtype=str)

for item in filepaths:
    convolve_data(item, resolution, resonumber)


# multiple lines

# here I create a duplicate of the HARPSfilelist and create the convolved HARPSfilelist  by replacing the ending S1D to S1Dresonumber to all its items


name_of_input = 'HARPSfilelist.dat'
name_of_output = 'res'+resonumber+name_of_input

input_list = open(name_of_input,'r')
output_list = open(name_of_output,'a')
for line in input_list:
    output_list.write(line.replace('S1D','S1Dres'+resonumber))
output_list.close()
input_list.close()

 
  
filepaths = np.loadtxt(name_of_output, dtype=str)

wavelength_range = np.loadtxt(linerange+'_lines.rdb', skiprows=2)
dw = 0.4
plot = False

directory_name = 'res'+resonumber+'_EWmyresults'

if not os.path.exists(directory_name):
    os.makedirs(directory_name)
    

try:
    size_of_filepaths = len(filepaths)
except TypeError:
    filepaths = [str(filepaths)]
    size_of_filepaths = len(filepaths)

    

for i in np.arange(size_of_filepaths):
    output = np.empty(len(wavelength_range))    
    for j in np.arange(len(wavelength_range)):
        output[j] = pseudo_EW(fname=filepaths[i], w1=wavelength_range[j,0], w2=wavelength_range[j,1], dw=dw, plot=plot)
    np.savetxt('./'+directory_name+'/result_'+filepaths[i].replace('.fits','.dat').replace('spectra/HARPS/',''), output, fmt='%.2f')
    

##### adding the 3 script


    
myData = "res"+resonumber+"_EWmyresults"
    

our_data, our_datanames = Load_Data_Pairs(myData)
    
    
wavelength_range = np.loadtxt(linerange+'_lines.rdb', skiprows=2)

wcentral = np.empty(len(wavelength_range))

for i in np.arange(len(wavelength_range)):
    winit = wavelength_range[i,0]
    wfin = wavelength_range[i,1]
    wcentral[i] = (winit + wfin)/2
np.savetxt(linerange+'_centrallines.dat', wcentral)
lines = np.loadtxt(linerange+"_centrallines.dat")
table = np.empty((len(lines), len(our_datanames)+1))
table[:,0] = lines
    
    
headers = "names"

for i in range(len(our_datanames)):
    headers = headers + "," + our_datanames[i]
    
    
    
for i in range(len(our_datanames)):
    table[:, 1+i] = our_data[i]
    print(len(our_data[i]))
        
np.savetxt("res"+resonumber+"_myEW.dat", table, header = headers, delimiter=",")   
    
     # transpose
     
table = np.loadtxt('res'+resonumber+'_myEW.dat', dtype=str, delimiter=',', comments="?")

    
table_T = np.transpose(table)
table_T[0,0] = table_T[0,0].replace('# ','')
print(table_T[0])
np.savetxt('res'+resonumber+'_transposed.csv', table_T, delimiter=',', fmt='%s')
print (table_T[:,0]) 
    
    # add the parameters to train and test
    
df = pd.read_csv("myoriginalparameters.dat", sep=" ", header = 0)

df.loc[:,"starname"] = df.Star.str.split("_").str[0].str.strip() 
#print (df.head())

df2 = pd.read_csv("res"+resonumber+"_transposed.csv", sep=",")

df2.loc[:,"starname"] = df2.names.str.split("_").str[1].str.strip()  # 1 assumes the files start with "results_"
#print (df2.head())
df3 = pd.merge(df2, df, on="starname", how="outer")
#df3[["names", "Star", "starname"]].head()
df3 = df3.drop(["starname", "Star"], axis=1)
#df3 = df2.join(df[["FeH", "Teff"]], on="starname")

    #print (df3.head())

df3.to_csv("res"+resonumber+"_linerangefor"+linerange+"_EWPar.csv", sep=",", index=False) 
