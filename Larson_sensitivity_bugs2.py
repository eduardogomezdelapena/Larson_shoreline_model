# -*- coding: utf-8 -*-
"""
Sensitivity analysis
"""

import matplotlib.pyplot as plt
import numpy as np
import math
from os.path import dirname, join as pjoin
import scipy.io as sio
import pandas as pd


#Initial conditions #Origin on inland

#Horizontal distances to origin
yl= 100 #Landward dune foot
dzid= 5 # for making the dune a trapezoid
ys= 120 #Seaward dune foot
yb= 170 #Berm crest seaward end

#Vertical distances
s= 5 #Dune crest height
Db= 1 #Berm crest height
Dc = 8 #Depth of closure

#Angles
beta= 3 # beach slope in degrees
###################################################################
yg= yb + math.tan((90-beta)*math.pi/180)*Db #Shoreline
#linear approach for subaqueous profile
yDc= yg + math.tan((90-beta)*math.pi/180)*Dc #Horizontal Dc
mwl= 0 #mean water level

####################################################################
#Load wave hindcast (Tairua)
hindfile= '/home/isabel/Escritorio/PhD_UoA/Shoreshop_data/Shoreshop_data/Case1/Waves/Wave_hindcast_corrected.mat'
mat_fname = pjoin(hindfile)

mat_contents = sio.loadmat(mat_fname)
sorted(mat_contents.keys())

mat_contents['hindcast'].dtype

time= mat_contents['hindcast'][0,0]['time']
H0= mat_contents['hindcast'][0,0]['Hs'] #Wave height
T=  mat_contents['hindcast'][0,0]['Tp'] #Wave period
L0= (9.8 * np.power(T,2)) / (2*math.pi) #Wave length, deep water

##############################################################################
#Define arrays

#Dune erosion and overwash coefficients
Cs= 0.001
A = 3

#Bar-berm material coefficients
Cb = 0.022
m= -0.5
lamb0= 0.15


#####################Transport equations##############################
######################Dune and erosion overwash#######################
# Cs= 1 * (10^-3)
#719529 is the number of days from matlab epoch to Unix epoch
#units are days
#humantime= pd.to_datetime(time[0]-719529,unit='d')
R= 0.158*np.sqrt(H0*L0) #Runup height
ZD= Db - mwl #distante from mwl to dune foot
#Cross shore transport rate

#The Runup height must exceed the dune foot
#elevation (R>ZD), otherwise there is no dune erosion.
whereqD=np.where((R>ZD)  == True)
print("Dune erosion (qD) occurs: "+str(np.sum(R>ZD))+" times")

qD= np.zeros(shape=(len(R),1))
for i in range(0,len(R)):
    if R[i]>ZD:
        qD[i]= 4 *Cs * ( (R[i] -ZD) *s / T[i]   )
    else:
        qD[i]=0

#Setting the overwash
#Where does R-ZD>s?
whereoverwash=np.where((R-ZD > s) == True)
print("Overwash occurs (qL): "+ str(np.sum(R-ZD > s))+" times")

alpha= np.zeros(shape=(len(R),1))
for i in range(0,len(R)):

    #For overwash to occur R-ZD > s
    if (R[i]-ZD) < s:
        alpha[i] = 0
    else:
        # A = 3
        alpha[i] = (1/A) * ( ((R[i]-ZD)/s) -1  )  


#Backwash
qs = qD / (1 + alpha)

#Overwash
qL= alpha*qD / (1+alpha)

######################Bar-berm material exchange#######################
# Cb = 1
# m= 1
# lamb0= 1
#Initial bar volume
VB0= 100 #m3 /m ???
w = 0.3 #sediment fall speed #Ferguson and Church 2004

#Equilibrium bar volume
VBE= np.power(L0,2) * Cb * (H0/(w*T)) * np.power(H0/L0,4/3)

lamb = lamb0 * np.power( (H0/(w*T)),m)

VB= VBE - (VBE-VB0)* np.exp(-lamb*(time-time[0]))

######################Wind-blown sand ################################
#Assuming a constant acretion rate of  2m3/yr
#Accretion rate= qwl + qws

qwl= 2 /365/24/8 / 2 #every 3 hours as in hindcast data
qws = qwl

######################################################################
#####################Conservation equations###########################

tend= 26813 #Number of iterations
dt = 0.125 #every 3 hours. units= day

qb= np.zeros(shape=(tend,1))

yb_arr= np.zeros(shape=(tend,1)); yb_arr[0]=yb
ys_arr= np.zeros(shape=(tend,1)); ys_arr[0]=ys
yl_arr= np.zeros(shape=(tend,1)); yl_arr[0]=yl
yg_arr= np.zeros(shape=(tend,1)); yg_arr[0]=yg


for i in range(1, tend): #len(time)):
    
    #Bar evolution
    qb[i] = (VB[i]-VB[i-1]) / dt 

    #Berm evolution
    yb_arr[i]= (dt * (-qws-qb[i] +qs[i])/(Db+Dc)) + yb_arr[i-1]
    
    #Dune evolution
    ys_arr[i]= ( dt  * (-qD[i] +qws)/s   ) + ys_arr[i-1]

    yl_arr[i] = ( dt  * (-qL[i] -qwl)/s   ) + yl_arr[i-1]
    #Shoreline evolution
    yg_arr[i] = yb_arr[i] + yg -yb 
#######################################################################

print("\n Parameters of interest:")
print("Cs:" + "%.6f" % Cs +", A:" + "%.2f" % A
      +", lambda0:"+ "%.2f" % lamb0 +", m:"+ "%.2f" % m +", Cb:"+ "%.2f" % Cb)
print("\n Transport rates:")
print("qD:" + str(qD[tend-1]) +", qL:" + str(qL[tend-1])
      +", qB:"+  str(qb[tend-1]))
print("\n State variables:")
print("yl:" + "%.2f" % yl_arr[-1] +", ys:" + "%.2f" % ys_arr[-1]
      +", yg:"+ "%.2f" % yg_arr[-1])
###################################################################
#State variables plot

plt.figure(1)
plt.plot(yg_arr[0:tend], color= "darkslategrey")
bottom,top= plt.ylim(); plt.ylim(bottom,top)
plt.xlabel("Time step")
plt.ylabel(" yg [meters]")
plt.title("Shoreline evolution " r'$(y_{g})$')

plt.figure(2)
plt.plot(ys_arr[0:tend], color= "indigo")
bottom,top= plt.ylim(); plt.ylim(bottom,top)
plt.xlabel("Time step")
plt.ylabel("ys [meters]")
plt.title("Dune seaward end " r'$(y_{s})$')

plt.figure(3)
plt.plot(yl_arr[0:tend], color= "orangered")
bottom,top= plt.ylim(); plt.ylim(bottom,top)
plt.xlabel("Time step")
plt.ylabel("yl [meters]")
plt.title("Dune landward end " r'$(y_{L})$')


#Transport rates plot
print("Dune erosion (qD) occurs: "+str(np.sum(qD[0:tend]>0))+" times")

plt.figure(4)
plt.plot(qD[0:tend], color= "mediumblue")
bottom,top= plt.ylim(); plt.ylim(bottom,top)
plt.xlabel("Time step")
plt.ylabel(r'$q_{D} [m^{3} s^{-1} m^{-1}]$')
plt.title("Cross-shore transport " r'$(q_{D})$')

print("Overwash (qL) occurs: "+str(np.sum(qL>0))+" times")
plt.figure(5)
plt.plot(qL, color= "indigo")
bottom,top= plt.ylim(); plt.ylim(bottom,top)
plt.xlabel("Time step")
plt.ylabel(r'$q_{L} [m^{3} s^{-1} m^{-1}]$')
plt.title("Overwash transport " r'$(q_{L})$')


plt.figure(6)
plt.plot(qb[0:tend], color= "orangered")
bottom,top= plt.ylim(); plt.ylim(bottom,top)
plt.xlabel("Time step")
plt.ylabel(r'$q_{b} [m^{3} s^{-1} m^{-1}]$')
plt.title("Bar-berm transport " r'$(q_{b})$')



