#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 14:48:54 2018

@author: aantoniadis
"""

from build_pew_verygeneral import pseudo_EW, read_data
import numpy as np
import os
from astropy.io import fits
from PyAstronomy import pyasl
from scipy.interpolate import interp1d

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


def convolve_data(fname, resolution):


    flux = fits.getdata(fname)
    hdr = fits.getheader(fname)
    w0, dw, N = hdr['CRVAL1'], hdr['CDELT1'], hdr['NAXIS1']
    wavelength = w0 + dw * np.arange(N)
    
    if round(dw,4) != 0.010:
        cdelt1 = 0.010
        f2 = interp1d(wavelength, flux, kind='linear')
        wavelength = np.arange(wavelength[0], wavelength[-1], cdelt1)
        flux = f2(wavelength)
        hdr['CDELT1']=0.010
    
    newflux = pyasl.instrBroadGaussFast(wavelength, flux, resolution, edgeHandling="firstlast", fullout=False, maxsig=None)
    
    fits.writeto(fname.replace('.fits', '')+'final'+'.fits', newflux, hdr, overwrite=True)


def EWmeasurements(convolution = "yes"):
    
    convo_limit= 100000
    spectralist = np.loadtxt('1Dfilelist.dat', dtype=str)
    filepaths = spectralist[:,0]
    resolution = spectralist[:,1]
    
    if convolution == "yes" :
    
        for i in np.arange(len(filepaths)):
            if np.float(resolution[i]) > convo_limit and np.float(resolution[i]) < 115000:
                convolve_data(filepaths[i], np.float(resolution[i]))
            else:
                os.system("cp "+filepaths[i]+" "+filepaths[i].replace('.fits', '')+'final'+'.fits')
        
        
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
        
        
          
        filepaths = np.loadtxt(name_of_output, dtype=str)
                
        name_of_input = '1Dfilelist.dat'
        name_of_output = 'final'+name_of_input
       
        input_list = open(name_of_input,'r')
        output_list = open(name_of_output,'w')
        
        for line in input_list:   
            output_list.write(line.replace('.fits','final'+'.fits'))
        output_list.close()
        
        input_list.close()
       
       
         
        filepaths = np.loadtxt(name_of_output, dtype=str)
        
        
        
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
                    
        
        for i in np.arange(size_of_filepaths):
            wavelength_range = np.loadtxt('lines.rdb', skiprows=2)
            
            print (wavelength_range)
            starwavelength, ___ =read_data(filepaths[i][0])
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
           
            # Load data and return the pairs
            our_data, our_datanames = Load_Data_Pairs(myData)
           
            table = np.empty((len(lines), 2))
            table[:,0] = lines
            table[:,1] = output
       
            if np.shape(filepaths[i])==(2,):
                  headers = "newstars" + "," + filepaths[i,0].replace('.fits','').replace('spectra/'+'newstars/','')
                  np.savetxt('./'+directory_name+filepaths[i,0].replace('.fits','').replace('spectra/'+'newstars/','')+"_EW.dat", table, header = headers, delimiter=",")
                  table = np.loadtxt('./'+directory_name+filepaths[i,0].replace('.fits','').replace('spectra/'+'newstars/','')+'_EW.dat', dtype=str, delimiter=',', comments="?")
            else:
                  headers = "newstars" + "," + filepaths[i,0].replace('.fits','').replace('spectra/'+'newstars/','')
                  np.savetxt('./'+directory_name+filepaths[i].replace('.fits','').replace('spectra/'+'newstars/','')+"_EW.dat", table, header = headers, delimiter=",")
                  table = np.loadtxt('./'+directory_name+filepaths[i].replace('.fits','').replace('spectra/'+'newstars/','')+'_EW.dat', dtype=str, delimiter=',', comments="?") 
      
       
           
            table_T = np.transpose(table)
            
            table_T[0,0] = table_T[0,0].replace('# ','')
            
            #print(table_T[0])
           
            if np.shape(filepaths[i])==(2,):
                  np.savetxt('./'+directory_name+filepaths[i,0].replace('.fits','').replace('spectra/'+'newstars/','')+'_newstars.csv', table_T, delimiter=',', fmt='%s')
            else:
                  np.savetxt('./'+directory_name+filepaths[i].replace('.fits','').replace('spectra/'+'newstars/','')+'_newstars.csv', table_T, delimiter=',', fmt='%s')
       
           
            #print (table_T[:,0])    
 

