#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 25 15:48:00 2018

@author: Matt
"""

from scipy import *
import numpy as np
from scipy.special import *
from time import time
from matplotlib.pyplot import *
import gaussxw as gs

# =============================================================================
# First we Will imprt functions that we will be useing
# calculate integrals using the 3 methods
# Plot the results
# =============================================================================


#define function to be integrated for erf
def f(x):
    """takes an x value and returns the 
    associated output of the function that
    gets integrated in the erf function"""
    return (2/sqrt(pi))*numpy.exp(-(x**2))

#define a function that integrates and function using trabezoidal rule
def trap(g,x1,x2,n):
    """takes a function g a start value of
    x1, as stop value of x2, and a number of 
    slices n and integrates useing the 
    trapezoidal method"""
    H=(x2-x1)/n#steps for integral
    v=0.5*(g(x1)+g(x2))#end of the trapezoidal rule steps aren counted as heavily so add seperate
    for j in range (1,n):#add all of the double counted internal trapezoids
        v+=(g(x1+j*H))
    return (v)*H

#define a function that integrates any function using simpsons rule
def simp(g,x1,x2,n):
# =============================================================================
#      """takes a function g a start value of
#     x1, as stop value of x2, and a number of 
#     slices n and integrates useing simpson's
#     rule"""
# =============================================================================
    h=(x2-x1)/n
    v=g(x1)+g(x2)#endpoint weighted differently add first
    for k in range(1,n):#i even
        if k%2==0:
            v+=2*g(x1+k*h)
        else:#i odd
            v+=4*g(x1+k*h)
    I=(1/3)*h*(v)#sum and multiply by 
    return I

def gaussint(g,x1,x2,n):
    x,w=gausswx(n)
    xp=0.5*(x2-x1)*x+0.5(x2+x1)
    wp=0.5*(x2-x1)
    ret=g(xp)
    ret=ret*wp
    ret=sum(ret)
    return ret

N=np.linspace(8,1000,1000-8+1,dtype=int)
print(N)
a=0
b=3
Itrap=np.empty(len(N))
Isimp=np.empty(len(N))
Igaus=np.empty(len(N))
for i in range (0,len(N)):
    Itrap[i]=trap(f,a,b,N[i])
    Isimp[i]=simp(f,a,b,N[i])
    Igaus[i]=gaussint(f,a,b,N[i])

