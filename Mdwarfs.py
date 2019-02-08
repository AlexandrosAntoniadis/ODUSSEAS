#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 15:15:35 2018

@author: aantoniadis
"""
import New_dataset 
import MachineLearning



######## SETTINGS ##############

option = "measure" # choose "measure" for both measuring the EWs of new stars and get parameters with ML, or "justML" to get results with different run/model 
regression = 'ridge' # choose the ML model: linear, ridge, ridgecv. Recommended: ridge
convo_limit = 100000 # if the resolution of the unknown spectra is higher than this, they will be convolved to the same resolution too.

########################


if option == "measure" :
    New_dataset.EWmeasurements(convo_limit)     
    
    MachineLearning.ML(regression)
    
elif option =="justML" :
    MachineLearning.ML(regression)
    
