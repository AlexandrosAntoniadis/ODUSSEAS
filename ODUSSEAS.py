#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 15:15:35 2018

@author: aantoniadis
"""
import New_dataset 
import New_dataset_withrv
import MachineLearning


######## SETTINGS ##############

option = "measure" # choose "measure" to measure the EWs of the new stars and get parameters with ML, or "justML" to get results with different run/model 
regression = 'ridge' # choose the ML model: linear, ridge, ridgecv, multitasklasso, multitaskelasticnet  (Recommended: ridge)
convo_limit = 100000 # if the resolution of the new spectra is higher than this, they will be convolved to their own resolution too.

########################


if option == "measure" :
    New_dataset.EWmeasurements(convo_limit)     
    
    MachineLearning.ML(regression)

elif option =="measure_withrv":
    New_dataset_withrv.EWmeasurements(convo_limit)     
    
    MachineLearning.ML(regression)

    
elif option =="justML" :
    MachineLearning.ML(regression)
    
