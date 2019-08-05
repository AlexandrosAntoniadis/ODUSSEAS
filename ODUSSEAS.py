#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 15:15:35 2018

@author: aantoniadis
"""

import New_data_withrv
import MachineLearning


######## SETTINGS ##############

regression = 'ridge' # choose the ML model: linear, ridge, ridgecv, multitasklasso, multitaskelasticnet  (Recommended: ridge)

########################

New_data_withrv.EWmeasurements()     
    
MachineLearning.ML(regression)

    

    
