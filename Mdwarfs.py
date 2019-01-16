#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 15:15:35 2018

@author: aantoniadis
"""
import New_dataset 
import MachineLearning

######### SETTINGS ##############

regression = 'ridge' # choose the ML model: linear, ridge, ridgecv. Recommended: ridge
inst = 'UVES' # put the name of the instrument of unknown stars
resolution = 110000 #set the resolution of unknown spectra
resonumber = '110' # set resolution number in k units
convolution_of_new = "yes" #convolve the unknown stars to their own resolution if it is comparable to the 115k of the HARPS dataset
linerange= "53to68" # the range of the linelist created for the spectra of the unknown stars. for example its total name should be: 53to68lines.rdb)
option = "justML" # choose "measure" for both measuring the EWs of new stars and get parameters with ML, or "justML" to get results with different run/model 

########################

if option == "measure" :
    New_dataset.EWmeasurements(convolution_of_new, resolution, resonumber, inst, linerange)     
    
    MachineLearning.ML(regression,convolution_of_new, resolution, resonumber, inst, linerange)
    
elif option =="justML" :
    MachineLearning.ML(regression,convolution_of_new, resolution, resonumber, inst, linerange)
    
    