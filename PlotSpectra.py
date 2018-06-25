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
          'figure.figsize': (10, 6),   #defines the default image size
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
plt.rcParams.update(params)



object_name1='ats_exp_0.1_20180615231436'
input_file1=object_name1+'_wl.txt'
figlog1='spec'+'_log_'+object_name1+'.png'
figlin1='spec'+'_lin_'+object_name1+'.png'

object_name2='ats_exp_0.5_20180615020051'
input_file2=object_name2+'_wl.txt'
figlog2='spec'+'_log_'+object_name2+'.png'
figlin2='spec'+'_lin_'+object_name2+'.png'


# identification by eye
pixel_order1=[482.89,518.741,651.33,688.739,912,967,1096.7]
lambda_order1=[404.656,435.883,546.074,578.013,763.51,811.53,912.3]
    
pixel_orderm1=[480.63,519.835,651.811,690.62,911.76,970.]
lambda_orderm1=[404.656,435.883,546.074,578.013,763.51,811.53]

# identification by eye
pixel=[-1200.,-970.0,-911.76,-690.62,-651.811,-519.835,-480.63,0,482.89,518.741,651.33,688.739,912,967,1096.7,1200.]
lambdawl=[-1000.,-811.53,-763.51,-578.013,-546.074,-435.883,-404.656,0,404.656,435.883,546.074,578.013,763.51,811.53,912.3,1000.]
    

#-------------------------------------------------
def plotspectrum(thefile,theobject,thelogfigname,thelinfigname):

    data=np.loadtxt(thefile)

    pixel=data[:,0]
    spec=data[:,1]


    plt.figure()
    plt.semilogy(pixel,spec,'b-')
    plt.grid(True)
    plt.xlabel('$\lambda$ (nm)',fontsize=15,fontweight='bold')
    plt.ylabel('ADU sum',fontsize=15,fontweight='bold')
    title='spectrum for {}'.format(theobject)
    plt.title(title,fontsize=20,fontweight='bold')
    plt.xlim(-1000.,1000.)
    
    plt.savefig(thelogfigname)
    plt.show()
    
    
    
    plt.figure()
    plt.plot(pixel,spec,'b-')
    plt.grid(True)
    plt.xlabel('$\lambda$ (nm)',fontsize=15,fontweight='bold')
    plt.ylabel('ADU sum',fontsize=15,fontweight='bold')
    title='spectrum for {}'.format(theobject)
    plt.title(title,fontsize=20,fontweight='bold')
    plt.xlim(-1000.,1000.)
    
    plt.savefig(thelinfigname)
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
    title='Dispersion relation for {}'.format(object_name1)
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
     title='Fit of Dispersion relation for {}'.format(object_name1)
     plt.title(title,fontsize=20,fontweight='bold')
     plt.show()
     
def Fitline2():
    
     x= pixel
     y=lambdawl 
     z = np.polyfit(x, y, 1)
     p = np.poly1d(z)
     
   
     
     theX=np.linspace(-1200.,1200,50)
     plt.figure()
     plt.plot(theX,p(theX),'r-')
   
     
     plt.scatter(pixel,lambdawl,marker='o',color='red')
     
     plt.grid(True)
     #plt.legend(loc='best')
     plt.xlabel('pixel',fontsize=15,fontweight='bold')
     plt.ylabel('$\lambda$ (nm)',fontsize=15,fontweight='bold')
     title='Fit of Dispersion relation for {}'.format(object_name1)
     plt.title(title,fontsize=20,fontweight='bold')
     plt.show()
    
    
    
    
    
if __name__ == "__main__":
    # execute only if run as a script
    plotspectrum(input_file1,object_name1,figlog1,figlin1)
    plotspectrum(input_file2,object_name2,figlog2,figlin2)
    calibrate()
    Fitline()
    Fitline2()

