#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 11 19:41:29 2018

@author: Matt
"""

#1) define constants and initial conditions given
#2) define functions for newtons law timesteps
#3) set up empty arrays for postion and velocity and add start points
#4) graph results for position and velocity

#import modules
from numpy import *
from scipy import *
from matplotlib import *
from pylab import *

#define constants and initial condition
MS=1#mass of the sun in kg
AU=1.496e11#m
G=6.67e-11#m**3/(kg*s**2)
GAU=39.5#AU**3/(MS*yr**2)
a=0.01#AU**2
x0=0.47#AU
y0=0#AU
vx0=0#au/yr
vy0=8.17#au/yr
dt=0.0001#timestep in years
T=1#total time
n=T/dt+1
n=int(n)

#define empty sets for 
X=zeros(n)
Y=zeros(n)
VX=zeros(n)
VY=zeros(n)
TT=linspace(0,T,n)
#set first index of arrays to initial conditions
X[0]=x0
Y[0]=y0
VX[0]=vx0
VY[0]=vy0

#define functions for newtons equations
#next nelocity
def Vnext(vj,xj,yj):
    r=((xj**2)+(yj**2))**0.5
    out=vj-((GAU*MS*xj)/((r)**3))*dt
    return out
#next position
def Xnext(vj,xj):
    out=xj+vj*dt
    return out

#run the loop
for i in range (0,n-1):
    VX[i+1]=Vnext(VX[i],X[i],Y[i])
    VY[i+1]=Vnext(VY[i],Y[i],X[i])
    X[i+1]=Xnext(VX[i+1],X[i])
    Y[i+1]=Xnext(VY[i+1],Y[i])
VV=(VX**2+VY**2)**0.5
#plot
figure(figsize=(6,10))
subplot(2,1,1)
plot(X,Y,label='orbit')
axis('equal')
xlabel('x position')
ylabel('y position')
title('orbit positon')
legend()
subplot(2,1,2)
plot(TT,VX,label='x velocity')
plot(TT,VY,label='y velocity')
plot(TT,VV,label='net velocity')
xlabel('time')
ylabel('velocity')
title('velocity')
legend()
tight_layout()