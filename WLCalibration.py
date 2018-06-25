#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 22 12:59:10 2018

@author: dagoret
"""


import numpy as np
import matplotlib
import matplotlib.pyplot as plt




# to enlarge the sizes
params = {'legend.fontsize': 'x-large',
          'figure.figsize': (15, 10),   #defines the default image size
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
plt.rcParams.update(params)



object_name='ats_exp_0.1_20180615231436'
input_file=object_name+'.txt'

# identification by eye
pixel_order1=[482.89,518.741,651.33,688.739,912,967,1096.7]
lambda_order1=[404.656,435.883,546.074,578.013,763.51,811.53,912.3]
    
pixel_orderm1=[480.63,519.835,651.811,690.62,911.76,970.]
lambda_orderm1=[404.656,435.883,546.074,578.013,763.51,811.53]




# identification by eye
pixel=[-1200.,-970.0,-911.76,-690.62,-651.811,-519.835,-480.63,482.89,518.741,651.33,688.739,912,967,1096.7,1200.]
lambdawl=[-1000.,811.53,-763.51,-578.013,-546.074,-435.883,-404.656,404.656,435.883,546.074,578.013,763.51,811.53,912.3,1000.]
    

#-------------------------------------------------
def plotspectrum():

    data=np.loadtxt(input_file)

    pixel=data[:,0]
    spec=data[:,1]


    plt.figure()
    plt.semilogy(pixel,spec,'b-')
    plt.grid(True)
    plt.xlabel('pixel',fontsize=15,fontweight='bold')
    plt.ylabel('ADU sum',fontsize=15,fontweight='bold')
    title='spectrum for {}'.format(object_name)
    plt.title(title,fontsize=20,fontweight='bold')
    plt.show()
    
    
    
    plt.figure()
    plt.plot(pixel,spec,'b-')
    plt.grid(True)
    plt.xlabel('pixel',fontsize=15,fontweight='bold')
    plt.ylabel('ADU sum',fontsize=15,fontweight='bold')
    title='spectrum for {}'.format(object_name)
    plt.title(title,fontsize=20,fontweight='bold')
    plt.show()
#-------------------------------------------------
def calibrate():
    

    
    plt.figure()
    plt.plot(pixel_order1,lambda_order1,'r:')
    plt.scatter(pixel_order1,lambda_order1,marker='o',color='red',label='order +1')
    
    plt.plot(pixel_orderm1,lambda_orderm1,'b:')
    plt.scatter(pixel_orderm1,lambda_orderm1,marker='o',color='blue',label='order -1')
    
    plt.grid(True)
    plt.legend(loc='best')
    plt.xlabel('pixel',fontsize=15,fontweight='bold')
    plt.ylabel('$\lambda$ (nm)',fontsize=15,fontweight='bold')
    title='Dispersion relation for {}'.format(object_name)
    plt.title(title,fontsize=20,fontweight='bold')
    plt.show()
#-------------------------------------------------------------
def Fitline():
    
     x= pixel_order1
     y=lambda_order1  
     z1 = np.polyfit(x, y, 1)
     p1 = np.poly1d(z1)
     
     x= pixel_orderm1
     y=lambda_orderm1  
     z2 = np.polyfit(x, y, 1)
     p2 = np.poly1d(z2)
     
     theX=np.linspace(0.,1200,50)
     plt.figure()
     plt.plot(theX,p1(theX),'r-')
     plt.plot(theX,p2(theX),'b-')
     
     plt.scatter(pixel_order1,lambda_order1,marker='o',color='red',label='order +1')
     plt.scatter(pixel_orderm1,lambda_orderm1,marker='o',color='blue',label='order -1')
     
     plt.grid(True)
     plt.legend(loc='best')
     plt.xlabel('pixel',fontsize=15,fontweight='bold')
     plt.ylabel('$\lambda$ (nm)',fontsize=15,fontweight='bold')
     title='Fit of Dispersion relation for {}'.format(object_name)
     plt.title(title,fontsize=20,fontweight='bold')
     plt.show()
     
    
if __name__ == "__main__":
    # execute only if run as a script
    plotspectrum()
    calibrate()
    Fitline()

