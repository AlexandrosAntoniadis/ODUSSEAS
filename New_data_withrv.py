#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 19:12:26 2019

@author: aantoniadis
"""

from build_pew_verygeneral import pseudo_EW, read_data
import numpy as np
import matplotlib.pyplot as plt
import os
from astropy.io import fits
from PyAstronomy import pyasl
from scipy.interpolate import interp1d

#####################

def gaussian(x, mu, sig):
  return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def find_rv(wavelength, flux, mask_lines=[5371.50, 5394.64, 5397.10, 5405.80, 5409.80, 5429.70, 5432.55, 5434.52, 5446.90, 5535.50, 6090.20, 6102.72, 6119.52, 6122.22], delta_rv=200):
  wave_mask = np.arange(mask_lines[0]-5., mask_lines[-1]+5., 0.01)
  flux_mask = 1. + wave_mask * 0.0
  for line in mask_lines:
    flux_mask -= gaussian(wave_mask, line, 0.1)

  ind_wi = np.where(wavelength < wave_mask[0])[0][-1]
  ind_wf = np.where(wavelength > wave_mask[-1])[0][0]

  wave_t = wavelength[ind_wi:ind_wf]
  flux_t = flux[ind_wi:ind_wf]


  rv, cc = pyasl.crosscorrRV(wave_t, flux_t, wave_mask, flux_mask, -delta_rv, delta_rv, 0.1, skipedge=500)
  maxind = np.argmax(cc)

 # plt.figure(1)
 # plt.plot(wave_t,flux_t/np.median(flux_t))
 # plt.plot(wave_mask,flux_mask)
 # plt.figure(2)
 # plt.plot(rv, cc)
 # plt.show()
 # plt.plot(rv, cc, 'bp-')
 # plt.plot(rv[maxind], cc[maxind], 'ro')
 # plt.show()
  
  return rv[maxind]

def rv_correction(fname, cdelt1=0.010):
    
    flux = fits.getdata(fname)
    hdr = fits.getheader(fname)
    w0, dw, N = hdr['CRVAL1'], hdr['CDELT1'], hdr['NAXIS1']
    wave = w0 + dw * np.arange(N)

    rv=find_rv(wave, flux)
    
    print ("RV:", rv)
    
    if abs(rv)>=1:
        
        wave_rv = wave / (1 + rv/3.e5)
        f2 = interp1d(wave_rv, flux, kind='linear')
        wave_int = np.arange(wave_rv[0], wave_rv[-1]-cdelt1, cdelt1)
        hdr['CRVAL1'] = wave_rv[0]

        flux_int = f2(wave_int)
        flux = flux_int
    else:
        wave_rv = wave

        hdr['CRVAL1'] = wave_rv[0]
    fits.writeto(fname, flux, hdr, overwrite=True)
    

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


def EWmeasurements(do_rv_cor):
  
        spectralist = np.loadtxt('1Dfilelist.dat', dtype=str)

        if spectralist.ndim > 1:
            filepaths = spectralist[:,0]
            resolution = spectralist[:,1]
        else:
        
            filepaths = [spectralist[0]]
            resolution = [spectralist[1]]

    
        for i in np.arange(len(filepaths)):
  
            os.system("cp "+filepaths[i]+" "+filepaths[i].replace('.fits', 'final')+''+'.fits')
                
    
            if do_rv_cor == 'yes':

               rv_correction(filepaths[i].replace('.fits', 'final')+''+'.fits')
                

        name_of_input = '1Dfilelist.dat'
        name_of_output = 'final'+name_of_input
        
        input_list = open(name_of_input,'r')
        output_list = open(name_of_output,'w')
        
        for line in input_list:    
            print(type(line))
            if np.shape(line) == (1,2):
                   output_list.write(line[:,0].replace('.fits','final'+'.fits'))
            else:    
                   output_list.write(line.replace('.fits','final'+'.fits'))
        output_list.close()
            
        input_list.close()
               
        name_of_input = '1Dfilelist.dat'
        name_of_output = 'final'+name_of_input
       
        input_list = open(name_of_input,'r')
        output_list = open(name_of_output,'w')
        
        for line in input_list:   
            output_list.write(line.replace('.fits','final'+'.fits'))
        output_list.close()
        
        input_list.close()
              
         
        filepaths = np.loadtxt(name_of_output, dtype=str, unpack=True, usecols=(0,))               
        print(filepaths)
        dw = 0.4
        plot = False
        
        directory_name = 'EWmyresults'
        
        if not os.path.exists(directory_name):
            os.makedirs(directory_name)
                    
        try:
            size_of_filepaths = len(filepaths)
        except TypeError:
            filepaths = [str(filepaths)]
            size_of_filepaths = len(filepaths)
        print(filepaths)              
        
        for i in np.arange(size_of_filepaths):
            wavelength_range = np.loadtxt('lines.rdb', skiprows=2)
            
            print (wavelength_range)
            starwavelength, ___ =read_data(filepaths[i])
            starwavelength_min, starwavelength_max = np.min(starwavelength), np.max(starwavelength)
            w1_max, w2_min = np.max(wavelength_range[:,0]) , np.min(wavelength_range[:,1]) 
            
            if starwavelength_min >  w2_min:
                map=np.where(wavelength_range[:,1] > starwavelength_min) 
                wavelength_range= wavelength_range[map[0]]
                
            if starwavelength_max <  w1_max:
                map=np.where(wavelength_range[:,0] < starwavelength_max) 
                wavelength_range= wavelength_range[map[0]]
            output = np.empty(len(wavelength_range))
            
        
            for j in np.arange(len(wavelength_range)):
                output[j] = pseudo_EW(fname=filepaths[i], w1=wavelength_range[j,0], w2=wavelength_range[j,1], dw=dw, plot=plot) 
                
            print(filepaths[i])
            
            if np.shape(filepaths[i])==(2,):
                np.savetxt('./'+directory_name+'/result_'+filepaths[i,0].replace('.fits','.dat').replace('spectra/'+'newstars/',''), output, fmt='%.2f')
              
            else:
                np.savetxt('./'+directory_name+'/result_'+filepaths[i].replace('.fits','.dat').replace('spectra/'+'newstars/',''), output, fmt='%.2f')
            
            wcentral = np.empty(len(wavelength_range))
        
            for k in np.arange(len(wavelength_range)):
                winit = wavelength_range[k,0]
                wfin = wavelength_range[k,1]
                wcentral[k] = (winit + wfin)/2
            
            if np.shape(filepaths[i])==(2,):
                  np.savetxt('./'+directory_name+filepaths[i,0].replace('.fits','').replace('spectra/'+'newstars/','')+'centrallines.dat', wcentral) 
                  lines = np.loadtxt('./'+directory_name+filepaths[i,0].replace('.fits','').replace('spectra/'+'newstars/','')+'centrallines.dat')
            else:
                  np.savetxt('./'+directory_name+filepaths[i].replace('.fits','').replace('spectra/'+'newstars/','')+'centrallines.dat', wcentral)  
                  lines = np.loadtxt('./'+directory_name+filepaths[i].replace('.fits','').replace('spectra/'+'newstars/','')+'centrallines.dat')
            

            myData = "EWmyresults" 
           
            our_data, our_datanames = Load_Data_Pairs(myData)
           
            table = np.empty((len(lines), 2))
            table[:,0] = lines
            table[:,1] = output
       
            if np.shape(filepaths[i])==(2,):
                  headers = "newstars" + "," + filepaths[i].replace('.fits','').replace('spectra/'+'newstars/','')
                  np.savetxt('./'+directory_name+filepaths[i].replace('.fits','').replace('spectra/'+'newstars/','')+"_EW.dat", table, header = headers, delimiter=",")
                  table = np.loadtxt('./'+directory_name+filepaths[i].replace('.fits','').replace('spectra/'+'newstars/','')+'_EW.dat', dtype=str, delimiter=',', comments="?")
            else:
                  headers = "newstars" + "," + filepaths[i].replace('.fits','').replace('spectra/'+'newstars/','')
                  np.savetxt('./'+directory_name+filepaths[i].replace('.fits','').replace('spectra/'+'newstars/','')+"_EW.dat", table, header = headers, delimiter=",")
                  table = np.loadtxt('./'+directory_name+filepaths[i].replace('.fits','').replace('spectra/'+'newstars/','')+'_EW.dat', dtype=str, delimiter=',', comments="?") 
      
       
            table_T = np.transpose(table)
            
            table_T[0,0] = table_T[0,0].replace('# ','')
            
           
            if np.shape(filepaths[i])==(2,):
                  np.savetxt('./'+directory_name+filepaths[i].replace('.fits','').replace('spectra/'+'newstars/','')+'_newstars.csv', table_T, delimiter=',', fmt='%s')
            else:
                  np.savetxt('./'+directory_name+filepaths[i].replace('.fits','').replace('spectra/'+'newstars/','')+'_newstars.csv', table_T, delimiter=',', fmt='%s')
       
