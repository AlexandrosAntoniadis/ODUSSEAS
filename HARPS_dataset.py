#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 31 12:54:57 2018

@author: aantoniadis
"""

import numpy as np
from astropy.io import fits
from PyAstronomy import pyasl
from build_pew import pseudo_EW

import os
import pandas as pd


######### SETTINGS ############## 

convolve_now = "yes" # convolve now if you want to add a new resolution dataset
resolution = 75000 # set the resolution that you want to convolve the HARPS stars.

#######################

resonumber = str(resolution)

c = 299792
DU_harps = c / 115000

R_new = resolution

DU_new = c / R_new

DU_conv = ((DU_new)**2 - (DU_harps)**2)**(0.5)

R_conv = c / (DU_conv)

def convolve_data(fname, R_conv, resonumber):


    flux = fits.getdata(fname)
    hdr = fits.getheader(fname)
    w0, dw, N = hdr['CRVAL1'], hdr['CDELT1'], hdr['NAXIS1']
    wavelength = w0 + dw * np.arange(N)
    
    newflux = pyasl.instrBroadGaussFast(wavelength, flux, R_conv, edgeHandling="firstlast", fullout=False, maxsig=None)
    
    fits.writeto(fname.replace('.fits', '')+'conv'+resonumber+'.fits', newflux, hdr, overwrite=True)

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



filepaths = np.loadtxt('RefHARPSfilelist.dat', dtype=str)

if convolve_now == "yes" :
    for item in filepaths:
        convolve_data(item, R_conv, resonumber)

     # below I create the convolved HARPSfilelist 

    name_of_input = 'RefHARPSfilelist.dat'
    name_of_output = 'res'+resonumber+name_of_input
    
    input_list = open(name_of_input,'r')
    output_list = open(name_of_output,'w')
    for line in input_list:
        output_list.write(line.replace('S1D.fits','S1Dconv'+resonumber+'.fits'))
    output_list.close()
    input_list.close()

else:
    name_of_output = 'RefHARPSfilelist.dat'

  
filepaths = np.loadtxt(name_of_output, dtype=str)

wavelength_range = np.loadtxt('lines.rdb', skiprows=2)
dw = 0.4
plot = False

directory_name = 'res'+resonumber+'RefHARPS_EWmyresults'

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
    

   
myData = "res"+resonumber+"RefHARPS_EWmyresults"
    
our_data, our_datanames = Load_Data_Pairs(myData)
    
    
wavelength_range = np.loadtxt('lines.rdb', skiprows=2)

wcentral = np.empty(len(wavelength_range))

for i in np.arange(len(wavelength_range)):
    winit = wavelength_range[i,0]
    wfin = wavelength_range[i,1]
    wcentral[i] = (winit + wfin)/2
np.savetxt('centrallines.dat', wcentral)
lines = np.loadtxt("centrallines.dat")
table = np.empty((len(lines), len(our_datanames)+1))
table[:,0] = lines
    
    
headers = "names"

for i in range(len(our_datanames)):
    headers = headers + "," + our_datanames[i]
       
    
for i in range(len(our_datanames)):
    table[:, 1+i] = our_data[i]
    
        
np.savetxt("res"+resonumber+"RefHARPS_myEW.dat", table, header = headers, delimiter=",")   
    

# transpose
     
table = np.loadtxt('res'+resonumber+'RefHARPS_myEW.dat', dtype=str, delimiter=',', comments="?")    
table_T = np.transpose(table)
table_T[0,0] = table_T[0,0].replace('# ','')
print(table_T[0])
np.savetxt('res'+resonumber+'RefHARPS_transposed.csv', table_T, delimiter=',', fmt='%s')
print (table_T[:,0]) 
    
# add the parameters as references for the ML to train and test
    
df = pd.read_csv("Refparameters.dat", sep=" ", header = 0)

df.loc[:,"starname"] = df.Star.str.split("_").str[0].str.strip() 


df2 = pd.read_csv("res"+resonumber+"RefHARPS_transposed.csv", sep=",")

df2.loc[:,"starname"] = df2.names.str.split("_").str[1].str.strip()  # since these files start with "result_"


df3 = pd.merge(df2, df, on="starname", how="outer")

df3 = df3.drop(["starname", "Star"], axis=1)

df3.to_csv("res"+resonumber+"_RefEWPar.csv", sep=",", index=False) 
