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
    
    
    newflux = pyasl.instrBroadGaussFast(wavelength, flux, resolution, edgeHandling="firstlast", fullout=False, maxsig=None)
    
    fits.writeto(fname.replace('.fits', '')+'final'+'.fits', newflux, hdr, overwrite=True)


def EWmeasurements(convo_limit, convolution_of_new = "yes"):

    spectralist = np.loadtxt('1Dfilelist.dat', dtype=str)
    filepaths = spectralist[:,0]
    resolution = spectralist[:,1]
    
    if convolution_of_new == "yes" :
    
        for i in np.arange(len(filepaths)):
            if np.float(resolution[i]) > convo_limit and np.float(resolution[i]) < 115000:
                convolve_data(filepaths[i], np.float(resolution[i]))
            else:
                os.system("cp "+filepaths[i]+" "+filepaths[i].replace('.fits', '')+'final'+'.fits')
        
        
        name_of_input = '1Dfilelist.dat'
        name_of_output = 'final'+name_of_input
        
        input_list = open(name_of_input,'r')
        output_list = open(name_of_output,'w')
        #for line in np.arange(len(filepaths)):
        #for line in input_list: #activate also lines 79 80
        for line in input_list:    
            print(type(line))
            if np.shape(line) == (1,2):
                   output_list.write(line[:,0].replace('.fits','final'+'.fits'))
            else:    
                   output_list.write(line.replace('.fits','final'+'.fits'))
        output_list.close()
        
        #    output_list.write(line.replace('.fits','final'+'.fits'))
        #output_list.close()
        #input_list.close()
        input_list.close()
        
        
          
        filepaths = np.loadtxt(name_of_output, dtype=str)
                
        name_of_input = '1Dfilelist.dat'
        name_of_output = 'final'+name_of_input
       
        input_list = open(name_of_input,'r')
        output_list = open(name_of_output,'w')
        #for line in np.arange(len(filepaths)):
        for line in input_list:   
            output_list.write(line.replace('.fits','final'+'.fits'))
        output_list.close()
        #input_list.close()
        input_list.close()
       
       
         
        filepaths = np.loadtxt(name_of_output, dtype=str)
        
        
        


        
        dw = 0.4
        plot = False
        
        directory_name = 'EWmyresults' #probably newstars_EWmyresults
        
        if not os.path.exists(directory_name):
            os.makedirs(directory_name)
                    
        try:
            size_of_filepaths = len(filepaths)
        except TypeError:
            filepaths = [str(filepaths)]
            size_of_filepaths = len(filepaths)
                    
        
        for i in np.arange(size_of_filepaths):
            wavelength_range = np.loadtxt('lines.rdb', skiprows=2)
            # indstead do : wavelength_rangeleft, wavelength_rangeright = np.genfromtxt('lines.rdb', names=['1','2']) #in this way it reads the two columns as seperate lists
            print (wavelength_range)
            starwavelength, ___ =read_data(filepaths[i][0])
            starwavelength_min, starwavelength_max = np.min(starwavelength), np.max(starwavelength)
            w1_max, w2_min = np.max(wavelength_range[:,0]) , np.min(wavelength_range[:,1]) 
            # instead do : w1_max, w2_min = np.max(wavelength_rangeleft) , np.min(wavelength_rangeright)
            if starwavelength_min >  w2_min:
                map=np.where(wavelength_range[:,1] > starwavelength_min) # instead of wavelength_range[j,1] . do it wavelength_rangeright
                wavelength_range= wavelength_range[map[0]]
                
            if starwavelength_max <  w1_max:
                map=np.where(wavelength_range[:,0] < starwavelength_max) # instead of wavelength_range[j,0] . do it wavelength_rangeleft
                wavelength_range= wavelength_range[map[0]]
            output = np.empty(len(wavelength_range))
            # the linelist is mapped but with no name or saved yet
        
            for j in np.arange(len(wavelength_range)):
                output[j] = pseudo_EW(fname=filepaths[i], w1=wavelength_range[j,0], w2=wavelength_range[j,1], dw=dw, plot=plot) 
                # above accordingly w1=wavelength_rangeleft and w2=wavelength_rangeleft
            print(filepaths[i])
            #np.savetxt('./'+directory_name+'/result_'+filepaths[i].replace('.fits','.dat').replace('spectra/'+'newstars/',''), output, fmt='%.2f') #when folder names change
            if np.shape(filepaths[i])==(2,):
                np.savetxt('./'+directory_name+'/result_'+filepaths[i,0].replace('.fits','.dat').replace('spectra/'+'newstars/',''), output, fmt='%.2f')
              
            else:
                np.savetxt('./'+directory_name+'/result_'+filepaths[i].replace('.fits','.dat').replace('spectra/'+'newstars/',''), output, fmt='%.2f')
            
            wcentral = np.empty(len(wavelength_range))
        
            for k in np.arange(len(wavelength_range)):
                winit = wavelength_range[k,0]
                wfin = wavelength_range[k,1]
                wcentral[k] = (winit + wfin)/2
            #np.savetxt('./'+directory_name+filepaths[i].replace('.fits','').replace('spectra/'+'newstars/','')+'centrallines.dat', wcentral)
            if np.shape(filepaths[i])==(2,):
                  np.savetxt('./'+directory_name+filepaths[i,0].replace('.fits','').replace('spectra/'+'newstars/','')+'centrallines.dat', wcentral) 
                  lines = np.loadtxt('./'+directory_name+filepaths[i,0].replace('.fits','').replace('spectra/'+'newstars/','')+'centrallines.dat')
            else:
                  np.savetxt('./'+directory_name+filepaths[i].replace('.fits','').replace('spectra/'+'newstars/','')+'centrallines.dat', wcentral)  
                  lines = np.loadtxt('./'+directory_name+filepaths[i].replace('.fits','').replace('spectra/'+'newstars/','')+'centrallines.dat')
            

#            #lines = np.loadtxt('./'+directory_name+filepaths[i].replace('.fits','').replace('spectra/'+'newstars/','')+'centrallines.dat')
            
            myData = "EWmyresults" #change according to folder (no inst, linerange)
           
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
            
            print(table_T[0])
           
            if np.shape(filepaths[i])==(2,):
                  np.savetxt('./'+directory_name+filepaths[i,0].replace('.fits','').replace('spectra/'+'newstars/','')+'_newstars.csv', table_T, delimiter=',', fmt='%s')
            else:
                  np.savetxt('./'+directory_name+filepaths[i].replace('.fits','').replace('spectra/'+'newstars/','')+'_newstars.csv', table_T, delimiter=',', fmt='%s')
       
           
            print (table_T[:,0])    
#            myData = "EWmyresults" #change according to folder (no inst, linerange)
#            
#            # Load data and return the pairs
#            our_data, our_datanames = Load_Data_Pairs(myData)
#            #print((len(lines), len(our_data)))
#            #table = np.empty((len(lines), len(our_datanames)+1))
#            #table[:,0] = lines
#        
#            
#            
#            headers = "newstars"
#            headers = headers + "," + our_datanames[i]
#            table = np.empty((len(our_data[i]), 2))
#            table[:,0] = lines
#            table[:, 1] = our_data[1]
#            
#            #for i in range(len(our_datanames)):
#               # headers = headers + "," + our_datanames[i]
#            
#            
#            
#            #for i in range(len(our_datanames)):
#                #table[:, 1+i] = our_data[i]
#                #print(len(our_data[i]))
#            
#            #np.savetxt('./'+directory_name+filepaths[i].replace('.fits','').replace('spectra/'+'newstars/','')+'final'+"_EW.dat", table, header = headers, delimiter=",")
#            
#            if np.shape(filepaths[i])==(2,):
#                  np.savetxt('./'+directory_name+filepaths[i,0].replace('.fits','').replace('spectra/'+'newstars/','')+'final'+"_EW.dat", table, header = headers, delimiter=",")
#            else:
#                  np.savetxt('./'+directory_name+filepaths[i].replace('.fits','').replace('spectra/'+'newstars/','')+'final'+"_EW.dat", table, header = headers, delimiter=",")
#            
#              
#            
#             # transpose
#             
#            #table = np.loadtxt('./'+directory_name+filepaths[i].replace('.fits','').replace('spectra/'+'newstars/','')+'final'+'_EW.dat', dtype=str, delimiter=',', comments="?")
#            
#            if np.shape(filepaths[i])==(2,):       
#                  table = np.loadtxt('./'+directory_name+filepaths[i,0].replace('.fits','').replace('spectra/'+'newstars/','')+'final'+'_EW.dat', dtype=str, delimiter=',', comments="?")
#            else:
#                  table = np.loadtxt('./'+directory_name+filepaths[i].replace('.fits','').replace('spectra/'+'newstars/','')+'final'+'_EW.dat', dtype=str, delimiter=',', comments="?") 
#            
#            table_T = np.transpose(table)
#            table_T[0,0] = table_T[0,0].replace('# ','')
#            print(table_T[0])
#            
#            
#            #np.savetxt('./'+directory_name+filepaths[i].replace('.fits','').replace('spectra/'+'newstars/','')+'final'+'_newstars.csv', table_T, delimiter=',', fmt='%s')
#            
#            if np.shape(filepaths[i])==(2,):
#                  np.savetxt('./'+directory_name+filepaths[i,0].replace('.fits','').replace('spectra/'+'newstars/','')+'final'+'_newstars.csv', table_T, delimiter=',', fmt='%s')
#            else:
#                  np.savetxt('./'+directory_name+filepaths[i].replace('.fits','').replace('spectra/'+'newstars/','')+'final'+'_newstars.csv', table_T, delimiter=',', fmt='%s')
#        
#            
#            print (table_T[:,0])
    
    

        
def EWmeasurementsAll(convolution_of_new, resolution, resonumber, inst, linerange):

    filepaths = np.loadtxt(inst+'1Dfilelist.dat', dtype=str)
    
    if convolution_of_new == "yes" :
    
        for item in filepaths:
            convolve_data(item, resolution, resonumber)
        
        
        name_of_input = inst+'1Dfilelist.dat'
        name_of_output = 'conv'+resonumber+name_of_input
        
        input_list = open(name_of_input,'r')
        output_list = open(name_of_output,'w')
        for line in input_list:
            output_list.write(line.replace('.fits','conv'+resonumber+'.fits'))
        output_list.close()
        input_list.close()
        
          
        filepaths = np.loadtxt(name_of_output, dtype=str)
                
        
        wavelength_range = np.loadtxt(linerange+'lines.rdb', skiprows=2)
        dw = 0.4
        plot = False
        
        directory_name = 'conv'+resonumber+inst+linerange+'_EWmyresults'
        
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
            
                
        myData = 'conv'+resonumber+inst+linerange+"_EWmyresults"
            
            # Load data and return the pairs
        our_data, our_datanames = Load_Data_Pairs(myData)
            
        
        wcentral = np.empty(len(wavelength_range))
        
        for i in np.arange(len(wavelength_range)):
            winit = wavelength_range[i,0]
            wfin = wavelength_range[i,1]
            wcentral[i] = (winit + wfin)/2
        np.savetxt(linerange+'centrallines.dat', wcentral)
        lines = np.loadtxt(linerange+'centrallines.dat')
        table = np.empty((len(lines), len(our_datanames)+1))
        table[:,0] = lines
            
            
        headers = "newstars"
        for i in range(len(our_datanames)):
            headers = headers + "," + our_datanames[i]
            
            
            
        for i in range(len(our_datanames)):
            table[:, 1+i] = our_data[i]
            print(len(our_data[i]))
                
        np.savetxt('conv'+resonumber+inst+linerange+"_EW.dat", table, header = headers, delimiter=",")   
            
             # transpose
             
        table = np.loadtxt('conv'+resonumber+inst+linerange+'_EW.dat', dtype=str, delimiter=',', comments="?")
        
        table_T = np.transpose(table)
        table_T[0,0] = table_T[0,0].replace('# ','')
        print(table_T[0])
        np.savetxt('conv'+resonumber+inst+"_linerange"+linerange+'_newstars.csv', table_T, delimiter=',', fmt='%s')
        
            
        print (table_T[:,0])
    
    
    else :
        
        
        
        wavelength_range = np.loadtxt(linerange+'lines.rdb', skiprows=2)
        dw = 0.4
        plot = False
        
        directory_name = inst+linerange+'_EWmyresults'
        
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
            
        
        
        myData = inst+linerange+"_EWmyresults"
            
            # Load data and return the pairs
        our_data, our_datanames = Load_Data_Pairs(myData)
            
        #wavelength_range = np.loadtxt(linerange+'_lines.rdb', skiprows=2)
        wcentral = np.empty(len(wavelength_range))
        
        for i in np.arange(len(wavelength_range)):
            winit = wavelength_range[i,0]
            wfin = wavelength_range[i,1]
            wcentral[i] = (winit + wfin)/2
        np.savetxt(linerange+'centrallines.dat', wcentral)
        lines = np.loadtxt(linerange+'centrallines.dat')
        table = np.empty((len(lines), len(our_datanames)+1))
        table[:,0] = lines
            
            
        headers = "newstars"
        for i in range(len(our_datanames)):
            headers = headers + "," + our_datanames[i]
            
            
            
        for i in range(len(our_datanames)):
            table[:, 1+i] = our_data[i]
            print(len(our_data[i]))
                
        np.savetxt(inst+linerange+"_EW.dat", table, header = headers, delimiter=",")   
            
             # transpose
             
        table = np.loadtxt(inst+linerange+'_EW.dat', dtype=str, delimiter=',', comments="?")
        
        table_T = np.transpose(table)
        table_T[0,0] = table_T[0,0].replace('# ','')
        print(table_T[0])
        np.savetxt(inst+"_linerange"+linerange+'_newstars.csv', table_T, delimiter=',', fmt='%s')
        
            
        print (table_T[:,0])      
                    
