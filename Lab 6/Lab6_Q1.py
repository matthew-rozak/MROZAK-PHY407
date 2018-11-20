#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 15 22:28:51 2018

@author: Matt
"""

#imports
import numpy as np
import matplotlib.pyplot as plt

#initial constant
G=1
M=10
L=2
#functions for dv/dt in x and y direction as given in the textbook
def der2x(x,y):
    """
    Takes the valuses x and y and computes their time derivative in 
    the x direction.
    """
    r=((x**2)+(y**2))**0.5
    ret=(-G*M*x)/((r**2)*((r**2)+(L**2)/4)**0.5)
    return ret
def der2y(x,y):
    """
    same as above but in the y direction
    """
    r=((x**2)+(y**2))**0.5
    ret=(-G*M*y)/((r**2)*((r**2)+(L**2)/4)**0.5)
    return ret

def f(r,t):
    """
    function to amalgamate the 4 derivatives for the integration
    """
    X=r[2]
    Y=r[3]
    vx=r[0]
    vy=r[1]
    dvx=der2x(X,Y)
    dvy=der2y(X,Y)
    dx=vx
    dy=vy
    return np.array([dvx,dvy,dx,dy],float)
#define time steps
a=0
b=10
N=1000
h=(b-a)/N
#define empty arrays to append things into
tpoints=np.arange(a,b,h)
xpoints=[]
ypoints=[]
vxpoints=[]
vypoints=[]
#define initial conditions
x=1
y=0
vx=0
vy=1
#set initial conditions into an array
r=np.array([vx,vy,x,y],float)
#integration loop
for t in tpoints:
    vxpoints.append(r[0])
    vypoints.append(r[1])
    xpoints.append(r[2])
    ypoints.append(r[3])
    k1 = h*f(r,t)
    k2 = h*f(r+0.5*k1,t+0.5*h)
    k3 = h*f(r+0.5*k2,t+0.5*h)
    k4 = h*f(r+k3,t+h)
    r += (k1+2*k2+2*k3+k4)/6
#plot
plt.scatter(xpoints,ypoints,c=tpoints,marker='.',cmap='viridis')
cbar=plt.colorbar()
cbar.set_label('time')
plt.xlabel('x position')
plt.ylabel('y position')
plt.axis("equal")
plt.title('orbit of space junk')
plt.savefig('orbit.png')
