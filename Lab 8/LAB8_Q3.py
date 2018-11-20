#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 31 10:53:49 2018

@author: Matt
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#constants
eps = 1
delx = 0.02
delt = 0.005
tf = 2
nts = int(tf/delt)
Lx = 2*np.pi
Tf = 2
bet = eps * delt / delx
X = np.arange(0,Lx,delx)
U = np.zeros([len(X),nts])
for i in range(1,len(U)-1):
    U[i][0] = np.sin(X[i])
for T in range(0,nts-1):
    if T == 0:
        for i in range(1,len(X)-1):
            U[i][T+1] = U[i][T]-(bet/2)*((U[i+1][T]**2)-(U[i-1][T]**2))
    else:
        for i in range(1,len(X)-1):
            U[i][T+1] = U[i][T-1]-(bet/2)*((U[i+1][T]**2)-(U[i-1][T]**2))
            
lowlim = -2
uplim = 2

plt.ion()

fig = plt.figure()
ax = fig.add_subplot(111)
line1, = ax.plot(X, U[:,0], 'r-')#,xlim=[0,1],ylim=[lowlim,uplim]) 
ax.set_xlim([0,Lx])
ax.set_ylim([lowlim,uplim])

for i in range(0,nts):
    if i%2==0:
        #ax.xlim([0,1])
        ax.set_ylim([lowlim,uplim])
        line1.set_ydata(U[:,i])
        fig.canvas.draw()
        fig.canvas.flush_events()

plt.show()

plt.subplot(221)
plt.ylim([lowlim,uplim])
plt.plot(X,U[:,int(0/delt)])
plt.xlabel('position')
plt.ylabel('displacment')
plt.title('time = %f'%(0))

plt.subplot(222)
plt.ylim([lowlim,uplim])
plt.plot(X,U[:,int(0.5/delt)])
plt.xlabel('position')
plt.ylabel('displacment')
plt.title('time = %f'%(0.5))

plt.subplot(223)
plt.ylim([lowlim,uplim])
plt.plot(X,U[:,int(1/delt)])
plt.xlabel('position')
plt.ylabel('displacment')
plt.title('time = %f'%(1))

plt.subplot(224)
plt.ylim([lowlim,uplim])
plt.plot(X,U[:,int(1.5/delt)])
plt.xlabel('position')
plt.ylabel('displacment')
plt.title('time = %f'%(1.5))

plt.tight_layout()
plt.savefig("Burgers.png")