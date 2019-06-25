#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 15:15:35 2018

@author: aantoniadis
"""
import New_data 
import New_data_withrv
import MachineLearning


######## SETTINGS ##############

option = "measure" # choose "measure" to measure the EWs of the new stars and get parameters with ML, or "justML" to get results with different run/model 
regression = 'ridge' # choose the ML model: linear, ridge, ridgecv, multitasklasso, multitaskelasticnet  (Recommended: ridge)

########################


if option == "measure" :
    New_data.EWmeasurements()     
    
    MachineLearning.ML(regression)

elif option == "justML" :
    MachineLearning.ML(regression)
    
    
#elif option =="measure_withrv":
#    New_data_withrv.EWmeasurements()     
    
#    MachineLearning.ML(regression)

    

    
