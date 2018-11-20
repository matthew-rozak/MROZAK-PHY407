#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 15 18:03:04 2018

@author: Matt
"""
"""
In order to solve the ODE you need to first define your acceleration functions
in the x and y direction as a vector so that it can be used with the integration
in a compact form. From here we are going to define a function that uses 
the verlet method to solve the coupled ODE's for the particles simultaneously.
This will be done by going through the components for each particle going
through equations 7 through 11 of the lab manual. We will include a time array 
so that we can have the plot colourmapped to this so we can see how the 
particle moves as tome progresses. From here we are going to plug in our 
initial conditions and use the function to give us our position arrays for each
particle. We will then plot the results.
"""
#2b
#imports
import numpy as np
import matplotlib.pyplot as plt
#define accleration functions for ax and ay
def ax(x,y):
    """
    takes x and y and calculates the derivative in the x direction
    """
    r=x**2+y**2
    if r==0:
        #this is necessary as then r=0 the particlebeing compared is itself
        #so there would be no acceleration due to itsself.
        return 0
    else:
        retx=2*((12*x/((r)**7))-(6*x/((r)**4)))
        return retx
def ay(x,y):
    """
    same as above but in the y direction
    """
    r=x**2+y**2
    if r==0:
        return 0
    else:
        rety=2*((12*y/((r)**7))-(6*y/((r)**4)))
        return rety
def a(x,y):
    """
    Combines above functions into one function
    """
    return ax(x,y), ay(x,y)

#define a function that solves the odes using verlet method 
def verl(A,V,R,ST,SPS,STN):
    """
    Takes a function that outputs x and y accaleration, V velocitie vectors
    R position vectors, ST start, SPS step size, and STN number of steps and
    integrates using verlet method and outputs x position and y position
    matrixes where the colums of the same number correspond to the matching 
    column in the other matrix.
    """
    #set initial conditions and define matrixes to set values for outputs
    x=np.empty([STN,len(R)])
    y=np.empty([STN,len(R)])
    vx=np.empty([STN,len(V)])
    vy=np.empty([STN,len(V)])
    #time arrays
    T=np.empty(STN)
    T[0]=ST
    h=SPS
    for i in range (0,len(R)):
        x[0][i]=R[i][0]
        y[0][i]=R[i][1]
        vx[0][i]=V[i][0]
        vy[0][i]=V[i][1]
    #loop for verlets method
    for i in range(0,STN-1):
        T[i+1]=T[i]+h#keep track of time
        if i == 0:#initial step corresponds to eq7
            VXt05h=vx[i]
            VYt05h=vy[i]
            for j in range (0,len(x[i])):
                for k in range(0,len(x[i])):
                    VXt05h += 1/2*h*A(x[i][j]-x[i][k],y[i][j]-y[i][k])[0]
                    VYt05h += 1/2*h*A(x[i][j]-x[i][k],y[i][j]-y[i][k])[1]
        #eq8
        x[i+1]=x[i]+h*VXt05h
        y[i+1]=y[i]+h*VYt05h
        #eq9
        kx=np.zeros(len(x[0]))
        ky=np.zeros(len(y[0]))
        for j in range (0,len(x[i])):
            for k in range(0,len(x[i])):
                kx[j]+=h*A(x[i+1][j]-x[i+1][k],y[i+1][j]-y[i+1][k])[0]
                ky[j]+=h*A(x[i+1][j]-x[i+1][k],y[i+1][j]-y[i+1][k])[1]
        #eq10
        vx[i+1]=VXt05h+0.5*kx
        vy[i+1]=VYt05h+0.5*ky
        #eq11
        VXt05h=VXt05h+kx
        VYt05h=VYt05h+ky
    return x,y,T#return arrays
#plots based off initial conditions
x,y,t=verl(a,[[0,0],[0,0]],[[4,4],[5.2,4]],0,0.01,100)
for i in range(len(x[0])):
    plt.scatter(x[:,i],y[:,i],c=t,marker='.',cmap='viridis')
cbar=plt.colorbar()
cbar.set_label('time')
plt.xlabel('x position')
plt.ylabel('y position')
plt.title('$\\vec{r}_1=[4,4]$,$\\vec{r}_2=[5.2,4]$')
plt.savefig('i.png')
plt.show()
x,y,t=verl(a,[[0,0],[0,0]],[[4.5,4],[5.2,4]],0,0.01,100)
for i in range(len(x[0])):
    plt.scatter(x[:,i],y[:,i],c=t,marker='.',cmap='viridis')
cbar=plt.colorbar()
cbar.set_label('time')
plt.xlabel('x position')
plt.ylabel('y position')
plt.title('$\\vec{r}_1=[4.5,4]$,$\\vec{r}_2=[5.2,4]$')
plt.savefig('ii.png')
plt.show()
x,y,t=verl(a,[[0,0],[0,0]],[[2,3],[3.5,4.4]],0,0.01,100)
for i in range(len(x[0])):
    plt.scatter(x[:,i],y[:,i],c=t,marker='.',cmap='viridis')
cbar=plt.colorbar()
cbar.set_label('time')
plt.xlabel('x position')
plt.ylabel('y position')
plt.title('$\\vec{r}_1=[2,3]$,$\\vec{r}_2=[3.5,4.4]$')
plt.savefig('iii.png')
plt.show()
