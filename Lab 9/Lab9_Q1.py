#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 09:26:41 2018

@author: Matt
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import dcst as dcst

#constants
L=1
d=0.1
C=1
sig=0.3
h=10**-3
n=90000
v=100
nl=4000
nts = 300
tf = h*nts
#define psi for initial conditions
def psi(x):
    a = C*x*(L-x)/L**2
    b = np.exp(-((x-d)**2)/(2*sig**2))
    return a*b
#matricies for time progression
ps = np.zeros([n,nts])
ph = np.zeros([n,nts])
l = np.linspace(0,L,n)#positions along the string
for i in range(0,n):
    ps[i,0] = psi(l[i])#initial velocities
#Fourier transform initial state
pstil = dcst.dst(ps[:,0])
phtil = dcst.dst(ph[:,0])
#Times we want
TT = (np.array([0,2,4,6,12,100]))*0.001
#integers for w values
nn = np.arange(0,n,1)
#calculate w to go along with each fourier term
w = np.pi*v*nn/L
for T in TT:#calculation at each time and plot
    sol = pstil/w*np.sin(w*T)
    sol = dcst.idst(sol)
    plt.plot(l,sol,label = T)
#lable plot
plt.legend(title = 'time (seconds)')
plt.xlabel('position')
plt.ylabel('displacment')
plt.tight_layout()
plt.savefig('spectral.png')
