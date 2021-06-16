# -*- coding: utf-8 -*-
"""
Editor de Spyder

Este es un archivo temporal
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
###########Make 2darray with coordinates###########################
coord= np.array([
    [0,yl,yl+dzid,ys-dzid,ys,yb,yg,yDc], #x coordinates
    [Db,Db,Db+s,Db+s,Db,Db,0,-Dc]  #y coordinates
    ])

mwl= 0 #mean water level 
sea = np.array([
    [yg,yDc],
    [mwl,mwl]
    ])
###################################################################
plt.plot(coord[0],coord[1],color="gold",marker="o")
plt.plot(sea[0],sea[1],color="dodgerblue")
plt.text(250, 7.5, "initial setting ", fontsize=12)
# plt.xlim(0, 200)
plt.ylim(-9,10)
plt.xlabel("Horizontal scale (m)")
plt.ylabel("Vertical scale (m)")
plt.title("Imaginary beach")
plt.show()

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

#####################Transport equations##############################
######################################################################
######################Dune and erosion overwash#######################
Cs= 1 * (10**-3)
#719529 is the number of days from matlab epoch to Unix epoch
#units are days
#humantime= pd.to_datetime(time[0]-719529,unit='d')
R= 0.158*np.sqrt(H0*L0) #Runup height
ZD= Db - mwl #distante from mwl to dune foot
#Cross shore transport rate
qD= 4 *Cs * ( np.power(R -ZD,s) / T   )
A = 3
alpha = (1/A) * ( ((R-ZD)/s) -1  )
#Backwash
qs = qD / (1 + alpha)
#Overwash
qL= alpha*qD / (1+alpha)

######################Bar-berm material exchange#######################
Cb = 1
m= 1
lamb0= 1
#Initial bar volume
VB0= 100 #m3 /m ¿?
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

dt = 0.125 #every 3 hours. units= day

tend= 100

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

    temp= np.array([
        [0,yl_arr[i],yl_arr[i]+dzid,ys_arr[i]-dzid,
         ys_arr[i],yb_arr[i],yg_arr[i],yDc], #x coordinates
        [Db,Db,Db+s,Db+s,Db,Db,0,-Dc]  #y coordinates
        ])
###################################################################

    plt.plot(temp[0],temp[1],color="gold",marker="o")
    plt.plot(sea[0],sea[1],color="dodgerblue")
    plt.text(250, 7.5, "time: "+str(i), fontsize=12)
    # plt.xlim(0, 200)
    plt.ylim(-9,10)
    plt.xlabel("Horizontal scale (m)")
    plt.ylabel("Vertical scale (m)")
    plt.title("Imaginary beach")
    plt.show()

plt.figure(2)
plt.plot(yg_arr[0:tend])
plt.xlabel("Time step")
plt.ylabel(" yg [meters]")
plt.title("Shoreline evolution")

plt.figure(3)
plt.plot(ys_arr[0:tend], color= "indigo")
plt.xlabel("Time step")
plt.ylabel("ys [meters]")
plt.title("Dune seaward end")

plt.figure(4)
plt.plot(yl_arr[0:tend], color= "orangered")
plt.xlabel("Time step")
plt.ylabel("yl [meters]")
plt.title("Dune landward end")
