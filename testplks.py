#!/usr/bin/env python

from __future__ import print_function

import time

t0 = time.time()

import emcee
#import triangle
import numpy as np
import scipy as sp
import scipy.stats
from numpy.linalg import inv
import scipy.optimize as op
import os.path as opa
import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator
import scipy.interpolate
import sys
#import Cmodules_runPB
import os
os.chdir('/Users/jgleyzes/Documents/Projets/EFTofLSS/Python/ResumCLASS/')
import zbeft, zbetools





Outpath='/Users/jgleyzes/Documents/Projets/EFTofLSS/Python/ResumCLASS/output/'
zbeftpath='/Users/jgleyzes/Documents/Projets/EFTofLSS/Python/ResumCLASS/'
classpath='/Users/jgleyzes/Documents/Projets/EFTofLSS/class'
resumpath='/Users/jgleyzes/Documents/Projets/EFTofLSS/Python/ResumCLASS/resum_data/'    
cosmorefpath='/Users/jgleyzes/Documents/Projets/EFTofLSS/Python/ResumCLASS/cosmo_ref/'
keys=['b1','b2','b3','b4','b5','b6','b7','b8','b9','b10','DM','kren','UseCosmoRef','ImportResummationMatrix','ExportResummationMatrix',
    'CLASS_path',#'outpath','zbEFT_path','PathToFolderRD','PathToFolderCosmoRef',
    'EpsAbs_NoCosmoRef' ,'EpsRel_NoCosmoRef' 
    ,'EpsAbs_YesCosmoRef' ,'EpsRel_YesCosmoRef' ,'h']
    
valuesbs=[ 2.150806600885515,
 -5.17990292025685,
 0.9723311553678212,
 6.845631078414692,
 6.765384087625802,
 -26.756867562986244,
 4.265894123752706,
 -54.085000005820774,
 116.453683753668,
 -81.51132551976502
,False,0.001,'no','no','yes',
classpath,#Outpath,zbeftpath,resumpath,cosmorefpath,
0.01,1e-5,0.1,5e-2,0.7 ]
pars = dict(zip(keys, valuesbs))
#kPSdata,PSdata,PSerr=np.loadtxt("/Users/jgleyzes/Documents/Projets/EFTofLSS/Python/ResumCLASS/output/ps1D_mean.dat").T
kPS=np.array([ 0.0160832,  0.0257545,  0.0356934,  0.045512 ,  0.0552959,
        0.0652349,  0.0751562,  0.0852165,  0.0952608,  0.105196 ,
        0.115136 ,  0.125083 ,  0.135091 ,  0.145084 ,  0.155103 ,
        0.165131 ,  0.175113 ,  0.185114 ,  0.195094 ,  0.205064 ,
        0.215101 ,  0.225118 ,  0.235075 ,  0.245071 ,  0.255076 ,
        0.265072 ,  0.275063 ,  0.285057 ,  0.295074 ,  0.305081 ,
        0.315064 ,  0.325061 ,  0.335065 ,  0.345049 ,  0.0160832,
       0.0257545,  0.0356934,  0.045512 ,  0.0552959,  0.0652349,
        0.0751562,  0.0852165,  0.0952608,  0.105196 ,  0.115136 ,
        0.125083 ,  0.135091 ,  0.145084 ,  0.155103 ,  0.165131 ,
        0.175113 ,  0.185114 ,  0.195094 ,  0.205064 ,  0.215101 ,
        0.225118 ,  0.235075 ,  0.245071 ,  0.255076 ,  0.265072 ,
        0.275063 ,  0.285057 ,  0.295074 ,  0.305081 ,  0.315064 ,
        0.325061 ,  0.335065 ,  0.345049 ,  0.0160832,  0.0257545,
        0.0356934,  0.045512 ,  0.0552959,  0.0652349,  0.0751562,
        0.0852165,  0.0952608,  0.105196 ,  0.115136 ,  0.125083 ,
        0.135091 ,  0.145084 ,  0.155103 ,  0.165131 ,  0.175113 ,
        0.185114 ,  0.195094 ,  0.205064 ,  0.215101 ,  0.225118 ,
        0.235075 ,  0.245071 ,  0.255076 ,  0.265072 ,  0.275063 ,
        0.285057 ,  0.295074 ,  0.305081 ,  0.315064 ,  0.325061 ,
        0.335065 ,  0.345049 ])



res=zbeft.plks_zbEFT(pars, kPS)
#np.savetxt('rescrefyesyesimph67.txt',res)
print(time.time()-t0)