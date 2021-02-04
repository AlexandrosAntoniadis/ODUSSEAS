#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 25 20:49:41 2018

@author: aantoniadis
"""

from __future__ import division
import numpy as np
import pandas as pd
import os
from time import time
import matplotlib.pyplot as plt
from sklearn import linear_model
from sklearn.model_selection import train_test_split
from sklearn.metrics import explained_variance_score, r2_score
#from sklearn.externals import joblib
import joblib

def mean_absolute_percentage_error(y_true, y_pred): 
    return np.mean(np.abs((y_true - y_pred) / y_true)) * 100

def mean_abso_error(y_true, y_pred):
    return np.mean(np.abs((y_true - y_pred)))

def mean_sq_error(y_true, y_pred):
    return np.mean(np.abs((y_true - y_pred)**2))
                                     

def ML(regression):
            
    if regression == 'linear':
        reg = linear_model.LinearRegression()
    elif regression == 'ridge':
        reg = linear_model.Ridge(alpha=1)
    elif regression == 'ridgecv': 
        reg = linear_model.RidgeCV(alphas=[0.1, 1.0, 10.0],cv=5,gcv_mode=None)
    elif regression == 'multitasklasso':
        reg = linear_model.MultiTaskLasso(alpha=1) 
    elif regression == 'multitaskelasticnet':
        reg = linear_model.MultiTaskElasticNet(alpha=1) 
        
             
    spectralist = np.loadtxt('final1Dfilelist.dat', dtype=str)
    if spectralist.ndim > 1:
        filepaths = spectralist[:,0]
        resolution = spectralist[:,1]
    else:
        filepaths = [spectralist[0]]
        resolution = [spectralist[1]]
    directory_name = 'EWmyresults'
    res_file=open('Parameter_Results.dat', 'w')
    res_file.write("# newstars [Fe/H] STD_FeH MAE_[Fe/H] Wide_Error_[Fe/H] Teff STD_Teff MAE_Teff Wide_Error_Teff R2_score EV_score \n")
                   
    
    MLplots_folder = 'Model_Prediction_Plots'

    if not os.path.exists(MLplots_folder):
        os.makedirs(MLplots_folder)  
             
    for i in np.arange(len(filepaths)):
                
        df = pd.read_csv('res'+resolution[i]+'_RefEWPar.csv')
        df.dropna(axis=1, inplace=True)
        
        names = ['names']
        
        y = df.columns.values[-2:]               
        
        newlines = np.loadtxt('./'+directory_name+filepaths[i].replace('.fits','').replace('spectra/'+'newstars/','')+'centrallines.dat',dtype=str)
        
        newdf = pd.read_csv('./'+directory_name+filepaths[i].replace('.fits','').replace('spectra/'+'newstars/','')+'_newstars.csv')
        
            
        newdf.dropna(axis=1, inplace=True)
        

        zeroews = newdf.mask(newdf<0.00001)
        newnames = ['newstars']
        newlabels = newdf[newnames]
        newdf = zeroews.dropna('columns')
        newlines = newdf.columns.values[1::]        

        
        df_x = df.drop(y,axis=1)
        
        df_y = df[y] 
                
                           
        df_x = pd.concat([df.loc[:, 'names'], df.loc[:, newlines]], axis=1)
                
        
        newdf = newdf.loc[:, newlines]
        
        FeH_list=[]
        Teff_list=[]
        MAE_FeH_list=[]
        MAE_Teff_list=[]
        Var_list=[]
        R2_list=[]
        
        for k in range(100): 
                  
            x_train, x_test, y_train, y_test = train_test_split(df_x, df_y, test_size=0.30)
            
            labels_train = x_train[names]
            labels_test = x_test[names]    
            
            x_train.drop(names,axis=1,inplace=True)
            x_test.drop(names,axis=1,inplace=True)
            
            reg.fit(x_train, y_train)
            
            joblib.dump(reg, 'savedmodel.pkl')
            
            y_pred_test = reg.predict(x_test)
            y_pred_train = reg.predict(x_train)
            
            
            N = len(y_test)
            starttime = time()
            y_pred_test = reg.predict(x_test)
            elapsedtime = time()
            t = elapsedtime - starttime
            
            
            #score_test = reg.score(x_test, y_test)
            
            variancescore = explained_variance_score(y_test, y_pred_test) 
            
            r2score = r2_score(y_test, y_pred_test)
            
            
            reg = joblib.load('savedmodel.pkl') #loading the saved model
            
            
            pred_newdf = reg.predict(newdf) #applying the saved model
            
            finalresults = pd.concat([newlabels, pd.DataFrame(pred_newdf)], axis=1)
            
            #print(finalresults)
            
            #print('Calculated parameters for {} stars in {:.2f}ms'.format(N, t*1e3))
            #print ('test score:', score_test)
            #print ('variance score:', variancescore)
            #print ('r2score:', r2score)
            
            mae_train = mean_abso_error(y_train[:],y_pred_train[:]) 
            mape_train = mean_absolute_percentage_error(y_train[:],y_pred_train[:]) 
            
            #print('Mean Absolute Error of Train : ' + str(mae_train))
            #print('Mean Absolute Percentage Error of Train : '+ str(mape_train))
            
            mae_test = mean_abso_error(y_test[:],y_pred_test[:])
            mape_test = mean_absolute_percentage_error(y_test[:],y_pred_test[:])
            
            #print('Mean Absolute Error of Test : ' + str(mae_test))
            #print('Mean Absolute Percentage Error of Test : '+ str(mape_test))
                    
            
            train_givenvalues = pd.concat([labels_train, y_train], axis=1)
            train_givenvalues = train_givenvalues.reset_index(drop=True)
            new_labeltrain=labels_train.reset_index(drop=True)
            train_predvalues = pd.concat([new_labeltrain, pd.DataFrame(y_pred_train)], axis=1)
                    
            test_givenvalues = pd.concat([labels_test, y_test], axis=1)
            test_givenvalues = test_givenvalues.reset_index(drop=True)
            new_labeltest=labels_test.reset_index(drop=True)
            test_predvalues = pd.concat([new_labeltest, pd.DataFrame(y_pred_test)], axis=1)
            
            
            FeH_list.append(finalresults[0])
            Teff_list.append(finalresults[1])
            MAE_FeH_list.append(mae_test[0])
            MAE_Teff_list.append(mae_test[1])
            R2_list.append(r2score)
            Var_list.append(variancescore)
            
        
        if resolution[i]=='115000':
            wide_error_FeH = ((np.std(FeH_list))**2 + (0.10)**2)**(1/2)
            wide_error_Teff = ((np.std(Teff_list))**2 + (65)**2)**(1/2)
        elif resolution[i]=='110000':
            wide_error_FeH = ((np.std(FeH_list))**2 + (0.10)**2)**(1/2)
            wide_error_Teff = ((np.std(Teff_list))**2 + (68)**2)**(1/2)
        elif resolution[i]=='94600':
            wide_error_FeH = ((np.std(FeH_list))**2 + (0.12)**2)**(1/2)
            wide_error_Teff = ((np.std(Teff_list))**2 + (77)**2)**(1/2)
        elif resolution[i]=='75000':
            wide_error_FeH = ((np.std(FeH_list))**2 + (0.13)**2)**(1/2)
            wide_error_Teff = ((np.std(Teff_list))**2 + (78)**2)**(1/2)
        elif resolution[i]=='48000':
            wide_error_FeH = ((np.std(FeH_list))**2 + (0.13)**2)**(1/2)
            wide_error_Teff = ((np.std(Teff_list))**2 + (80)**2)**(1/2)         
        
        
        res_file.write(str(newlabels.iat[0,0])+' '+ str(round(np.mean(FeH_list),3))+' '+ str(round(np.std(FeH_list),3))+' '+ str(round(np.mean(MAE_FeH_list),3))+' '+ str(round((wide_error_FeH),3))+' '+ str(int(np.mean(Teff_list)))+' '+ str(int(np.std(Teff_list)))+' '+ str(int(np.mean(MAE_Teff_list)))+' '+ str(int(wide_error_Teff))+' '+ str(round(np.mean(R2_list),3))+' '+ str(round(np.mean(Var_list),3)) + "\n")
        
        starname = filepaths[i].replace('.fits','').replace('spectra/'+'newstars/','')
  
        print('Star '+str(newlabels.iat[0,0])+' results completed and saved in Paremeter_Results.dat')
        
                   
        # plots of the FeH test
           
        set_res = 15
        fig, ax = plt.subplots(figsize=(set_res*0.8,set_res*0.5))
        ax.set_title('[Fe/H]'+' '+'model'+' '+'testing', fontsize=set_res*1.5)
        ax.set_xlabel("M.L. [Fe/H] [dex]", fontsize=set_res*1.5)
        ax.set_ylabel("Ref. [Fe/H] [dex]", fontsize=set_res*1.5)        
        ax.plot((-0.8,0.4),(-0.8,0.4),'--b', lw=2) # for FeH
        ax.plot(test_predvalues.values[:,1],test_givenvalues.values[:,1],'ko')      
        ax.tick_params(axis='both',labelsize=set_res*1.5)
        ax.spines['right'].set_visible(True)
        ax.spines['top'].set_visible(True)
        plt.savefig("./"+MLplots_folder+"/"+starname+'_FeH_test_comparison.pdf', bbox_inches='tight')
        plt.close()
        
        # Feh diff
        set_res = 15
        fig, ax = plt.subplots(figsize=(set_res*0.8,set_res*0.2))
        ax.set_xlabel("M.L. [Fe/H] [dex]", fontsize=set_res*1.5)
        ax.set_ylabel("$\Delta$[Fe/H] [dex]", fontsize=set_res*1.5)        
        ax.plot((-0.8,0.4),(0,0),'--b', lw=2) # for FeH 
        ax.plot(test_predvalues.values[:,1],test_predvalues.values[:,1]-test_givenvalues.values[:,1] ,'ko')      
        ax.tick_params(axis='both',labelsize=set_res*1.5)
        ax.spines['right'].set_visible(True)
        ax.spines['top'].set_visible(True)
        plt.savefig("./"+MLplots_folder+"/"+starname+'_Diff_FeH_test_comparison.pdf', bbox_inches='tight')
        plt.close()
        
        # plots of the Teff test
        set_res = 15
        fig, ax = plt.subplots(figsize=(set_res*0.8,set_res*0.5))
        ax.set_title('T$_{\mathrm{eff}}$'+' '+'model'+' '+'testing', fontsize=set_res*1.5)
        ax.set_xlabel("M.L. T$_{\mathrm{eff}}$ [K]", fontsize=set_res*1.5)
        ax.set_ylabel("Ref. T$_{\mathrm{eff}}$ [K]", fontsize=set_res*1.5)
        ax.tick_params(axis='both',labelsize=set_res*1.5)
        ax.spines['right'].set_visible(True)
        ax.spines['top'].set_visible(True)
        ax.plot((2700,4000),(2700,4000),'--b', lw=2) #for Teff
        ax.plot(test_predvalues.values[:,2], test_givenvalues.values[:,2],'ko')       
        ax.plot(clip_box=True, clip_on=True)
        plt.savefig("./"+MLplots_folder+"/"+starname+'_Teff_test_comparison.pdf', bbox_inches='tight')
        plt.close()
        
        # T diff
        set_res = 15
        fig, ax = plt.subplots(figsize=(set_res*0.8,set_res*0.2))
        ax.set_xlabel("M.L. T$_{\mathrm{eff}}$ [K]", fontsize=set_res*1.5)
        ax.set_ylabel("$\Delta$T$_{\mathrm{eff}}$ [K]", fontsize=set_res*1.5)
        ax.tick_params(axis='both',labelsize=set_res*1.5)
        ax.spines['right'].set_visible(True)
        ax.spines['top'].set_visible(True)
        ax.plot((2700,4000),(0,0),'--b', lw=2) #for Teff
        ax.plot(test_predvalues.values[:,2],test_predvalues.values[:,2]-test_givenvalues.values[:,2],'ko')       
        ax.plot(clip_box=True, clip_on=True)
        plt.savefig("./"+MLplots_folder+"/"+starname+'_Diff_Teff_test_comparison.pdf', bbox_inches='tight')
        plt.close()
           
        
        # plots of the FeH train
           
        #set_res = 15
        #fig, ax = plt.subplots(figsize=(set_res*0.8,set_res*0.5))
        #ax.set_title('[Fe/H]'+' '+'model'+' '+'training', fontsize=set_res*1.5)
        #ax.set_xlabel("M.L. [Fe/H] [dex]", fontsize=set_res*1.5)
        #ax.set_ylabel("Ref. [Fe/H] [dex]", fontsize=set_res*1.5)        
        #ax.plot((-0.8,0.4),(-0.8,0.4),'--b', lw=2) # for FeH
        #ax.plot(train_predvalues.values[:,1],train_givenvalues.values[:,1],'ko')      
        #ax.tick_params(axis='both',labelsize=set_res*1.5)
        #ax.spines['right'].set_visible(True)
        #ax.spines['top'].set_visible(True)
        #plt.savefig("./"+MLplots_folder+"/"+starname+'_FeH_train_comparison.pdf', bbox_inches='tight')
        #plt.close()
        
        
        # plots of the Teff train
        
        #set_res = 15
        #fig, ax = plt.subplots(figsize=(set_res*0.8,set_res*0.5))
        #ax.set_title('T$_{\mathrm{eff}}$'+' '+'model'+' '+'training', fontsize=set_res*1.5)
        #ax.set_xlabel("M.L. T$_{\mathrm{eff}}$ [K]", fontsize=set_res*1.5)
        #ax.set_ylabel("Ref. T$_{\mathrm{eff}}$ [K]", fontsize=set_res*1.5)
        #ax.tick_params(axis='both',labelsize=set_res*1.5)
        #ax.spines['right'].set_visible(True)
        #ax.spines['top'].set_visible(True)
        #ax.plot((2700,4000),(2700,4000),'--b', lw=2) #for Teff
        #ax.plot(train_predvalues.values[:,2],train_givenvalues.values[:,2],'ko')       
        #ax.plot(clip_box=True, clip_on=True)
        #plt.savefig("./"+MLplots_folder+"/"+starname+'_Teff_train_comparison.pdf', bbox_inches='tight')
        #plt.close()
        
    res_file.close()        
