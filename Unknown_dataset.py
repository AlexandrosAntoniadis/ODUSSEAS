#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 31 14:33:22 2018

@author: mac
"""

from build_pew_verygeneral import pseudo_EW
import numpy as np
import os
from astropy.io import fits
from PyAstronomy import pyasl

######################
# set the name of instrument. 
# set resolution and resonumber that need to be convolved
# set the range of the linelist you are working with

resolution = 110000
resonumber = "110"
inst = "UVES"
linerange = "53to58"

#####################
def Load_Data_Pairs(data1):
    
    # Get the data and the names from our results
    our_datanames = []
    our_data = []
    for file in os.listdir(data1):
        # get all .dat files
        if file.endswith(".dat"):
            temp_name = os.path.join(data1, file)
            # remove unwanded elements and append the names
            our_datanames.append(temp_name.split('/')[1].strip( '.dat' ))
            # append data
            our_data.append(np.loadtxt(temp_name))
    
        
    return our_data, our_datanames

def convolve_data(fname, resolution, resonumber):


    flux = fits.getdata(fname)
    hdr = fits.getheader(fname)
    w0, dw, N = hdr['CRVAL1'], hdr['CDELT1'], hdr['NAXIS1']
    wavelength = w0 + dw * np.arange(N)
    
    
    newflux = pyasl.instrBroadGaussFast(wavelength, flux, resolution, edgeHandling="firstlast", fullout=False, maxsig=None)
    
    #prihdr = fits.Header()
    #prihdr.set("CDELT1",w0)
    #prihdr.set("CRVAL1",dw)
    #prihdr.set("NAXIS1",N)

    fits.writeto(fname.replace('.fits', '')+'conv'+resonumber+'.fits', newflux, hdr, overwrite=True)



filepaths = np.loadtxt(inst+'1Dfilelist.dat', dtype=str)

for item in filepaths:
    convolve_data(item, resolution, resonumber)


name_of_input = inst+'1Dfilelist.dat'
name_of_output = 'conv'+resonumber+name_of_input

input_list = open(name_of_input,'r')
output_list = open(name_of_output,'a')
for line in input_list:
    output_list.write(line.replace('.fits','conv'+resonumber+'.fits'))
output_list.close()
input_list.close()

#name_of_output = 'conv'+resonumber+inst+'1Dfilelist.dat' 
#the above line should be setted manually, and commenting out the commands 72-80, in case of overcalculating due to multiple paths after consecutive runs
  
filepaths = np.loadtxt(name_of_output, dtype=str)

wavelength_range = np.loadtxt(linerange+'_lines.rdb', skiprows=2)
dw = 0.4
plot = False

directory_name = 'conv'+resonumber+inst+'_EWmyresults'

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
    np.savetxt('./'+directory_name+'/result_'+filepaths[i].replace('.fits','.dat').replace('spectra/'+inst+'1D/',''), output, fmt='%.2f')
    


myData = 'conv'+resonumber+inst+"_EWmyresults"
    
    # Load data and return the pairs
our_data, our_datanames = Load_Data_Pairs(myData)
    
wavelength_range = np.loadtxt(linerange+'_lines.rdb', skiprows=2)
wcentral = np.empty(len(wavelength_range))

for i in np.arange(len(wavelength_range)):
    winit = wavelength_range[i,0]
    wfin = wavelength_range[i,1]
    wcentral[i] = (winit + wfin)/2
np.savetxt(linerange+'_centrallines.dat', wcentral)
lines = np.loadtxt(linerange+'_centrallines.dat')
table = np.empty((len(lines), len(our_datanames)+1))
table[:,0] = lines
    
    
headers = "newstars"
for i in range(len(our_datanames)):
    headers = headers + "," + our_datanames[i]
    
    
    
for i in range(len(our_datanames)):
    table[:, 1+i] = our_data[i]
    print(len(our_data[i]))
        
np.savetxt('conv'+resonumber+inst+"_EW.dat", table, header = headers, delimiter=",")   
    
     # transpose
     
table = np.loadtxt('conv'+resonumber+inst+'_EW.dat', dtype=str, delimiter=',', comments="?")

table_T = np.transpose(table)
table_T[0,0] = table_T[0,0].replace('# ','')
print(table_T[0])
np.savetxt('conv'+resonumber+inst+'_newstars.csv', table_T, delimiter=',', fmt='%s')

    
print (table_T[:,0])
