#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 15:15:35 2018

@author: aantoniadis
"""

import New_data
import MachineLearning


######## SETTINGS ##############

regression = 'ridge' # choose the ML model: linear, ridge, ridgecv, multitasklasso, multitaskelasticnet  (Recommended: ridge)
do_rv_cor = 'yes' # choose 'yes' if you want to do rv correction to the spectra or 'no' if they are already corrected
########################

New_data.EWmeasurements(do_rv_cor)     
    
MachineLearning.ML(regression)
