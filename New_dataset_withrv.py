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

def find_rv(wavelength, flux, mask_lines=[6081.40, 6085.23, 6090.20, 6102.71, 6108.12, 6111.70, 6122.21, 6162.17], delta_rv=200):
  wave_mask = np.arange(mask_lines[0]-5., mask_lines[-1]+5., 0.01)
  flux_mask = 1. + wave_mask * 0.
  for line in mask_lines:
    flux_mask -= gaussian(wave_mask, line, 0.1)

  ind_wi = np.where(wavelength < wave_mask[0])[0][-1]
  ind_wf = np.where(wavelength > wave_mask[-1])[0][0]

  wave_t = wavelength[ind_wi:ind_wf]
  flux_t = flux[ind_wi:ind_wf]


  rv, cc = pyasl.crosscorrRV(wave_t, flux_t, wave_mask, flux_mask, -delta_rv, delta_rv, 0.1, skipedge=500)
  maxind = np.argmax(cc)

  plt.figure(1)
  plt.plot(wave_t,flux_t)
  plt.plot(wave_mask,flux_mask)
  plt.figure(2)
  plt.plot(rv, cc)
  plt.show()
#  plt.plot(rv, cc, 'bp-')
#  plt.plot(rv[maxind], cc[maxind], 'ro')
#  plt.show()
  return rv[maxind]

def rv_correction(fname, cdelt1=0.010):
    
    flux = fits.getdata(fname)
    hdr = fits.getheader(fname)
    w0, dw, N = hdr['CRVAL1'], hdr['CDELT1'], hdr['NAXIS1']
    wave = w0 + dw * np.arange(N)
    
    rv=find_rv(wave, flux)
    
    print ("RV:", rv)
    wave_rv = wave / (1 + rv/3.e5)
    f2 = interp1d(wave_rv, flux, kind='linear')
    wave_int = np.arange(wave_rv[0], wave_rv[-1], cdelt1)
    hdr['CRVAL1'] = wave_rv[0]
    flux_int = f2(wave_int)
    fits.writeto(fname, flux_int, hdr, overwrite=True)
    

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
    
    
    #plt.plot(wavelength, flux)
    #plt.show()
    
    newflux = pyasl.instrBroadGaussFast(wavelength, flux, resolution, edgeHandling="firstlast", fullout=False, maxsig=None)
    
    fits.writeto(fname.replace('.fits', '')+'final'+'.fits', newflux, hdr, overwrite=True)


def EWmeasurements(convo_limit, convolution = "yes",rv_cor=True):
    
    spectralist = np.loadtxt('1Dfilelist.dat', dtype=str)
#    print(np.shape(spectralist))
#    if np.shape(spectralist)==(2,):
#        filepaths = spectralist[0]
#        resolution = spectralist[1]
#    else:
#    print(type(spectralist))
#    print(len(spectralist))
#    input()
    filepaths = spectralist[:,0]
    resolution = spectralist[:,1]
    
    if convolution == "yes" :
    
        for i in np.arange(len(filepaths)):
#            print(i,len(filepaths),np.shape(filepaths))
            if np.float(resolution[i]) > convo_limit and np.float(resolution[i]) < 115000:
                convolve_data(filepaths[i], np.float(resolution[i]))
            else:
                os.system("cp "+filepaths[i]+" "+filepaths[i].replace('.fits', '')+'final'+'.fits')
    
            if rv_cor == True:
                rv_correction(filepaths[i].replace('.fits', '')+'final'+'.fits')
                #input()
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

    
    

        
#def EWmeasurementsAll(convolution_of_new, resolution, resonumber, inst, linerange):
#
#    filepaths = np.loadtxt(inst+'1Dfilelist.dat', dtype=str)
#    
#    if convolution_of_new == "yes" :
#    
#        for item in filepaths:
#            convolve_data(item, resolution, resonumber)
#        
#        
#        name_of_input = inst+'1Dfilelist.dat'
#        name_of_output = 'conv'+resonumber+name_of_input
#        
#        input_list = open(name_of_input,'r')
#        output_list = open(name_of_output,'w')
#        for line in input_list:
#            output_list.write(line.replace('.fits','conv'+resonumber+'.fits'))
#        output_list.close()
#        input_list.close()
#        
#          
#        filepaths = np.loadtxt(name_of_output, dtype=str)
#                
#        
#        wavelength_range = np.loadtxt(linerange+'lines.rdb', skiprows=2)
#        dw = 0.4
#        plot = False
#        
#        directory_name = 'conv'+resonumber+inst+linerange+'_EWmyresults'
#        
#        if not os.path.exists(directory_name):
#            os.makedirs(directory_name)
#                    
#        try:
#            size_of_filepaths = len(filepaths)
#        except TypeError:
#            filepaths = [str(filepaths)]
#            size_of_filepaths = len(filepaths)
#                    
#        
#        for i in np.arange(size_of_filepaths):
#            output = np.empty(len(wavelength_range))    
#            for j in np.arange(len(wavelength_range)):
#                output[j] = pseudo_EW(fname=filepaths[i], w1=wavelength_range[j,0], w2=wavelength_range[j,1], dw=dw, plot=plot)
#            np.savetxt('./'+directory_name+'/result_'+filepaths[i].replace('.fits','.dat').replace('spectra/'+inst+'1D/',''), output, fmt='%.2f')
#            
#                
#        myData = 'conv'+resonumber+inst+linerange+"_EWmyresults"
#            
#            # Load data and return the pairs
#        our_data, our_datanames = Load_Data_Pairs(myData)
#            
#        
#        wcentral = np.empty(len(wavelength_range))
#        
#        for i in np.arange(len(wavelength_range)):
#            winit = wavelength_range[i,0]
#            wfin = wavelength_range[i,1]
#            wcentral[i] = (winit + wfin)/2
#        np.savetxt(linerange+'centrallines.dat', wcentral)
#        lines = np.loadtxt(linerange+'centrallines.dat')
#        table = np.empty((len(lines), len(our_datanames)+1))
#        table[:,0] = lines
#            
#            
#        headers = "newstars"
#        for i in range(len(our_datanames)):
#            headers = headers + "," + our_datanames[i]
#            
#            
#            
#        for i in range(len(our_datanames)):
#            table[:, 1+i] = our_data[i]
#            print(len(our_data[i]))
#                
#        np.savetxt('conv'+resonumber+inst+linerange+"_EW.dat", table, header = headers, delimiter=",")   
#            
#             # transpose
#             
#        table = np.loadtxt('conv'+resonumber+inst+linerange+'_EW.dat', dtype=str, delimiter=',', comments="?")
#        
#        table_T = np.transpose(table)
#        table_T[0,0] = table_T[0,0].replace('# ','')
#        print(table_T[0])
#        np.savetxt('conv'+resonumber+inst+"_linerange"+linerange+'_newstars.csv', table_T, delimiter=',', fmt='%s')
#        
#            
#        print (table_T[:,0])
#    
#    
#    else :
#        
#        
#        
#        wavelength_range = np.loadtxt(linerange+'lines.rdb', skiprows=2)
#        dw = 0.4
#        plot = False
#        
#        directory_name = inst+linerange+'_EWmyresults'
#        
#        if not os.path.exists(directory_name):
#            os.makedirs(directory_name)
#            
#        
#        try:
#            size_of_filepaths = len(filepaths)
#        except TypeError:
#            filepaths = [str(filepaths)]
#            size_of_filepaths = len(filepaths)
#        
#            
#        
#        for i in np.arange(size_of_filepaths):
#            output = np.empty(len(wavelength_range))    
#            for j in np.arange(len(wavelength_range)):
#                output[j] = pseudo_EW(fname=filepaths[i], w1=wavelength_range[j,0], w2=wavelength_range[j,1], dw=dw, plot=plot)
#            np.savetxt('./'+directory_name+'/result_'+filepaths[i].replace('.fits','.dat').replace('spectra/'+inst+'1D/',''), output, fmt='%.2f')
#            
#        
#        
#        myData = inst+linerange+"_EWmyresults"
#            
#            # Load data and return the pairs
#        our_data, our_datanames = Load_Data_Pairs(myData)
#            
#        #wavelength_range = np.loadtxt(linerange+'_lines.rdb', skiprows=2)
#        wcentral = np.empty(len(wavelength_range))
#        
#        for i in np.arange(len(wavelength_range)):
#            winit = wavelength_range[i,0]
#            wfin = wavelength_range[i,1]
#            wcentral[i] = (winit + wfin)/2
#        np.savetxt(linerange+'centrallines.dat', wcentral)
#        lines = np.loadtxt(linerange+'centrallines.dat')
#        table = np.empty((len(lines), len(our_datanames)+1))
#        table[:,0] = lines
#            
#            
#        headers = "newstars"
#        for i in range(len(our_datanames)):
#            headers = headers + "," + our_datanames[i]
#            
#            
#            
#        for i in range(len(our_datanames)):
#            table[:, 1+i] = our_data[i]
#            print(len(our_data[i]))
#                
#        np.savetxt(inst+linerange+"_EW.dat", table, header = headers, delimiter=",")   
#            
#             # transpose
#             
#        table = np.loadtxt(inst+linerange+'_EW.dat', dtype=str, delimiter=',', comments="?")
#        
#        table_T = np.transpose(table)
#        table_T[0,0] = table_T[0,0].replace('# ','')
#        print(table_T[0])
#        np.savetxt(inst+"_linerange"+linerange+'_newstars.csv', table_T, delimiter=',', fmt='%s')
#        
#            
#        print (table_T[:,0])      
                    
