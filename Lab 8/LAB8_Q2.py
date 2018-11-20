#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 30 22:07:35 2018

@author: Matt
"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#%matplotlib qt

#constants
L=1
d=0.1
C=1
sig=0.3
h=10**-6
n=200
v=100
nl=400
nts = 30000
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
for T in range(0,nts-1):#integration loop
    for i in range(0,n):
        ph[i][T+1] = ph[i][T]+h*ps[i][T]#update positions
    for i in range(1,n-1):#hold edgepoints constant for velocities
        ps[i][T+1] = ps[i][T]+h*(v**2/(l[1]**2))*(ph[i+1][T]+ph[i-1][T]-2*ph[i][T])     
    if T%1000 == 0:
        print(T/nts*100,"%")#counter to see progressiosn
#lowlim = min(min(j) for j in ph)
#uplim = max(max(j) for j in ph)
lowlim = -0.001#define limits
uplim = 0.001#define limits
# =============================================================================
# fig, ax = plt.figure(), plt.axes(xlim=(0,1),ylim=(lowlim,uplim))
# x = l
# line, = ax.plot(x,ph[:,i])
# def init():  # only required for blitting to give a clean slate.
#     line.set_ydata([np.nan] * len(x))
#     return line,
# 
# 
# def animate(i):
#     line.set_ydata(ph[:,i])  # update the data.
#     return line,
# 
# 
# ani = animation.FuncAnimation(
#     fig, animate, init_func=init, interval=20, blit=False, frames = len(ph[3,:]))
# #ani.save("movie.mp4")
# 
# plt.show()
# =============================================================================
# =============================================================================
# xdata, ydata = [], []
# fig, ax = plt.figure(), plt.axes(xlim=(0,1),ylim=(lowlim,uplim))
# line, = ax.plot([], [], lw=2)
# def init():
#     line.set_data([], [])
#     return line,
# def animate(i):
#     x = l
#     y = ph[:,i]
#     line.set_data(x, y)
#     return line,
# 
# anim = FuncAnimation(fig, animate, frames=np.linspace(0, tf, nts/100),
#                     init_func=init, blit=True)
# anim.save('basic_animation.mp4',writer='ffmpeg',fps=50,dpi=600)
# =============================================================================
#make sure backend graphics are set to automatic
plt.ion()
#initial plot set up
fig = plt.figure()
ax = fig.add_subplot(111)
line1, = ax.plot(l, ph[:,0], 'r-')#,xlim=[0,1],ylim=[lowlim,uplim]) 
ax.set_xlim([0,1])
ax.set_ylim([lowlim,uplim])

for i in range(0,nts):
    if i%10==0:#skip some steps to make go faster.
        #ax.xlim([0,1])
        #updates plots at time progresses
        ax.set_ylim([lowlim,uplim])
        line1.set_ydata(ph[:,i])
        fig.canvas.draw()
        fig.canvas.flush_events()
plt.show()

plt.subplot(221)
plt.ylim([lowlim,uplim])
plt.plot(l,ph[:,3000])
plt.xlabel('position')
plt.ylabel('displacment')
plt.title('time = %f'%(h*3000))

plt.subplot(222)
plt.ylim([lowlim,uplim])
plt.plot(l,ph[:,15000])
plt.xlabel('position')
plt.ylabel('displacment')
plt.title('time = %f'%(h*15000))

plt.subplot(223)
plt.ylim([lowlim,uplim])
plt.plot(l,ph[:,23800])
plt.xlabel('position')
plt.ylabel('displacment')
plt.title('time = %f'%(h*23800))

plt.subplot(224)
plt.ylim([lowlim,uplim])
plt.plot(l,ph[:,28500])
plt.xlabel('position')
plt.ylabel('displacment')
plt.title('time = %f'%(h*28500))
plt.tight_layout()
plt.savefig('destabilization.png')
