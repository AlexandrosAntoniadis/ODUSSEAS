#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 25 20:49:41 2018

@author: mac
"""

from __future__ import division
import numpy as np
import pandas as pd
from time import time
import matplotlib.pyplot as plt
from sklearn import linear_model
from sklearn.model_selection import train_test_split
from sklearn.metrics import explained_variance_score, r2_score
from sklearn.externals import joblib

def mean_absolute_percentage_error(y_true, y_pred): 
    return np.mean(np.abs((y_true - y_pred) / y_true)) * 100

def mean_abso_error(y_true, y_pred):
    return np.mean(np.abs((y_true - y_pred)))

def mean_sq_error(y_true, y_pred):
    return np.mean(np.abs((y_true - y_pred)**2))
                                     

def ML(regression,convolution_of_new, resolution, resonumber, inst, linerange):
            
    newlines = np.loadtxt(linerange+'centrallines.dat',dtype=str)    
        
    df = pd.read_csv('conv'+resonumber+'_goodEWPar.csv')
    df.dropna(axis=1, inplace=True)
    
    
    if convolution_of_new == "yes" :
        newdf = pd.read_csv('conv'+resonumber+inst+"_linerange"+linerange+'_newstars.csv')
    else :
        newdf = pd.read_csv(inst+"_linerange"+linerange+'_newstars.csv')    
        
    newdf.dropna(axis=1, inplace=True) 
    
    
    if regression == 'linear':
        reg = linear_model.LinearRegression()
    elif regression == 'ridge':
        reg = linear_model.Ridge(alpha=1)
    elif regression == 'ridgecv': 
        reg = linear_model.RidgeCV(alphas=[0.1, 1.0, 10.0],cv=10,gcv_mode=None)
    elif regression == 'lassolars':
        reg = linear_model.LassoLars(alpha=1)
    
    
    names = ['names']
    
    y = df.columns.values[-2:] 
    labels = df[names]
    
    newnames = ['newstars']
    newlabels = newdf[newnames] 
    
    df_x = df.drop(y,axis=1)    
    
    df_x = pd.concat([df.loc[:, 'names'], df.loc[:, newlines]], axis=1)
            
    newdf = newdf.drop(newlabels,axis=1) 
    
    newdf = newdf.loc[:, newlines]
    
    
    df_y = df[y] 
        
    x_train, x_test, y_train, y_test = train_test_split(df_x, df_y, test_size=0.30)
    
    labels_train = x_train[names]
    labels_test = x_test[names]    
    
    x_train.drop(names,axis=1,inplace=True)
    x_test.drop(names,axis=1,inplace=True)
    
    reg.fit(x_train, y_train)
    
    #for saving the model run the next line. For keeping current,comment it out
    joblib.dump(reg, 'res'+resonumber+'savedmodel.pkl')
    
    y_pred_test = reg.predict(x_test)
    y_pred_train = reg.predict(x_train)
    
    
    N = len(y_test)
    starttime = time()
    y_pred_test = reg.predict(x_test)
    elapsedtime = time()
    t = elapsedtime - starttime
    
    
    score_test = reg.score(x_test, y_test)
    
    variancescore = explained_variance_score(y_test, y_pred_test) 
    
    r2score = r2_score(y_test, y_pred_test)
    
    
    reg = joblib.load('res'+resonumber+'savedmodel.pkl') #for loading the saved model
    
    newdf = reg.predict(newdf) #applying the saved model to the new stars
    
    
    print('Calculated parameters for {} stars in {:.2f}ms'.format(N, t*1e3))
    print ('test score:', score_test)
    print ('variance score:', variancescore)
    print ('r2score:', r2score)
    
    mae_train = mean_abso_error(y_train[:],y_pred_train[:]) 
    mape_train = mean_absolute_percentage_error(y_train[:],y_pred_train[:]) 
    
    print('Mean Absolute Error of Train : ' + str(mae_train))
    print('Mean Absolute Percentage Error of Train : '+ str(mape_train))
    
    mae_test = mean_abso_error(y_test[:],y_pred_test[:])
    mape_test = mean_absolute_percentage_error(y_test[:],y_pred_test[:])
    
    print('Mean Absolute Error of Test : ' + str(mae_test))
    print('Mean Absolute Percentage Error of Test : '+ str(mape_test))
    
    
    finalresults = pd.concat([newlabels, pd.DataFrame(newdf)], axis=1)
    
    print(finalresults)
    
    train_givenvalues = pd.concat([labels_train, y_train], axis=1)
    train_givenvalues = train_givenvalues.reset_index(drop=True)
    new_labeltrain=labels_train.reset_index(drop=True)
    train_predvalues = pd.concat([new_labeltrain, pd.DataFrame(y_pred_train)], axis=1)
    
    
    test_givenvalues = pd.concat([labels_test, y_test], axis=1)
    test_givenvalues = test_givenvalues.reset_index(drop=True)
    new_labeltest=labels_test.reset_index(drop=True)
    test_predvalues = pd.concat([new_labeltest, pd.DataFrame(y_pred_test)], axis=1)
    
    
    # plotting the train FeH
    set_res = 12
    
    plt.figure(figsize=([set_res,set_res]))
    plt.title('[Fe/H]'+' '+'train'+' '+'comparison')
    plt.ylabel("ML value", fontsize=set_res)
    plt.xlabel("AA value", fontsize=set_res)
    plt.plot((-0.8,0.8),(-0.8,0.8),'--b') # for Fe/H
    plt.plot(train_givenvalues.values[:,1], train_predvalues.values[:,1],'o')
    plt.show()    
    
    # plotting the train Teff
    set_res = 12
    
    plt.figure(figsize=([set_res,set_res]))
    plt.title('Teff'+' '+'train'+' '+'comparison')
    plt.ylabel("ML value", fontsize=set_res)
    plt.xlabel("AA value", fontsize=set_res)
    plt.plot((2000,4000),(2000,4000),'--b') #for Teff    
    plt.plot(train_givenvalues.values[:,2], train_predvalues.values[:,2],'o')    
    plt.show()
        
    # plotting the test FeH
    set_res = 12
    
    plt.figure(figsize=([set_res,set_res]))
    plt.title('[Fe/H]'+' '+'test'+' '+'comparison')
    plt.ylabel("ML value", fontsize=set_res)
    plt.xlabel("AA value", fontsize=set_res)
    plt.plot((-0.8,0.8),(-0.8,0.8),'--b') # for FeH
    plt.plot(test_givenvalues.values[:,1], test_predvalues.values[:,1],'o')
    plt.show()    
    
    # plotting the test Teff
    set_res = 12
    
    plt.figure(figsize=([set_res,set_res]))
    plt.title('Teff'+' '+'test'+' '+'comparison')
    plt.ylabel("ML value", fontsize=set_res)
    plt.xlabel("AA value", fontsize=set_res)
    plt.plot((2000,4000),(2000,4000),'--b') #for Teff
    plt.plot(test_givenvalues.values[:,2], test_predvalues.values[:,2],'o')
    plt.show()
    