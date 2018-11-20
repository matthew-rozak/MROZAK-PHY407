#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 22:32:19 2018

@author: Matt
"""
"""
Description of what the code beneth does:
In this code we are going to be using runge kutta 4 integration here in rk4
with constant stepsize and varied step size to integrate the equations of 
motion for a peice of space junk. PLease consider adaptive, varied and adjusted
to be synonyms. We will refer to rk4 with constant step size as rk4c, and rk4 
with varied step size as rk4v. There are also components to analyze the time it
takes each of the methods to run, and one to visualize how the timestep changes
throughout time in the adaptive stepsize RK4.
"""

import numpy as np
import matplotlib.pyplot as plt
from time import time

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

def rk4_step(h0,ra,t,g):
    """
    This takes a timestep, position and velocity vector, time, and 
    amalgamated derivative munction and outputs the ammount a rk4 
    function increases ra by
    """
    k1a = h0*g(ra,t)
    k2a = h0*g(ra+0.5*k1a,t+0.5*h0)
    k3a = h0*g(ra+0.5*k2a,t+0.5*h0)
    k4a = h0*g(ra+k3a,t+h0)
    ka = (k1a+2*k2a+2*k3a+k4a)/6
    return ka

#starting rk4c   
#define time steps for rk4c
a=0
b=10
N=10000
h=(b-a)/N
#define empty arrays to append things into for rk4c
tpoints=np.arange(a,b,h)
xpoints=[]
ypoints=[]
vxpoints=[]
vypoints=[]
#define initial conditions for both rk4v and rk4c
x=1
y=0
vx=0
vy=1
#set initial conditions into an array for rk4c
r=np.array([vx,vy,x,y],float)
#integration loop for rk4c
start=time()#start timer for rk4c loop
for t in tpoints:
    vxpoints.append(r[0])
    vypoints.append(r[1])
    xpoints.append(r[2])
    ypoints.append(r[3])
    r += rk4_step(h,r,t,f)
end=time()#end timmer for rk4c loop
diff=end-start#time taken for rk4c loop
print("time taken for integration useing Runge-Kutta method of integration with constant step size is",diff,"seconds")#prints time taken for rk4c
#plot for rk4c
plt.scatter(xpoints,ypoints,c=tpoints,marker='.',cmap='viridis')
cbar=plt.colorbar()
cbar.set_label('time rk4 constant time step')
plt.xlabel('x position')
plt.ylabel('y position')
plt.axis("equal")
plt.title('orbit of space junk')
#plt.savefig('orbit.png')
#set initial values of variables
#starting rk4v
#define initial conditions matrix for rk4v
ra=np.array([vx,vy,x,y],float)
t=[]#blank list to append numbers for rk4v
ti=0#initial time for rk4v
i=0#counter fo indicies of t  list for rk4v
delta=10**-6#error we want to be in rk4v
#returned lists from rk4v initlization
xapoints=[]
yapoints=[]
vxapoints=[]
vyapoints=[]
h=0.01#initial tiome step for rk4v
h0=h#set h0 to initial time step for rk4v
j=0#in the while loop set initial valuse into results matricies for rk4v
b=10#define time to go until in rk4v loop
#timer
start=time() #start timer for rk4v loop
#use a while loop for the adaptive time step as we dont know what the stepsize
#at each step will be, this allows us to go until a desired end time.
while ti < b:
    #use j==0 as a condition to append values from last run into results list
    #if j==0 then we have accepted our values as rho was greater than 1 if not
    #then we are reedoing the step with a smaller timestep in hopes of 
    #increasing rho. 
    if j==0:#accepts calculation from last ittteration when rho > 1
        t.append(ti)
        vxapoints.append(ra[0])
        vyapoints.append(ra[1])
        xapoints.append(ra[2])
        yapoints.append(ra[3])
    #first rk4 with h0
    raa = ra + rk4_step(h0,ra,t[i],f)
    #second rk4 with h0
    raa = raa + rk4_step(h0,raa,t[i]+h0,f)
    #rk4 with h1 = 2h0
    h1 = 2*h0
    rab = ra + rk4_step(h1,ra,t[i],f)
    #error between 2 rk4's with h0 and rk4 with h1=2*h0
    ex=(1/30)*(raa[2]-rab[2])
    ey=(1/30)*(raa[3]-rab[3])
    #calculate rho
    rho=h0*delta/((ex**2+ey**2)**0.5)
    #update h0
    h0 = h0*rho**(1/4)
    #condition where rho is greater than 1 and we accept or previous calculation
    if rho > 1:
        ra = raa#accpts value rk2 with 2 steps
        ti=ti+h1#increases time counter
        i+=1#incteases t index counter
        j=0#acccepts calculated values at beginning of next itteration
    else:
        j=1#rejects values and runs again with smaller h
end=time()# ends timer
diff=end-start#calculates time takes for rk4v
print("time taken for integration useing Runge-Kutta method of integration with variable step size is",diff,"seconds")#prints the time taken for rk4v
#add plot for rk4v 
plt.scatter(xapoints,yapoints,c=t,marker='.',cmap='magma_r') 
cbar=plt.colorbar()
cbar.set_label('time rk4 adaptive time step')
plt.xlabel('x position')
plt.ylabel('y position')
plt.axis("equal")
plt.title('orbit of space junk')
plt.savefig('adjustedtimesteporbit.png')
plt.show()
#calculating timestep used as time progresses
dtpoints=np.array(t[1:])-np.array(t[:-1])#calculates the difference between 
#consecutive times from the rk4v integration
plt.plot(t[:-1],dtpoints)
plt.xlabel('t time(s)')
plt.ylabel('$\\Delta t$ timestep(s)')
plt.title('how the timestep used changes as time progresses')
plt.savefig('timestepvariation.png')