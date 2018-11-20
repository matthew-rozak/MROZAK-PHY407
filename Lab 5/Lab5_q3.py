#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 20:52:12 2018

@author: Matt
"""
#a
#imports
import numpy as np
import matplotlib.pyplot as plt
import numpy.fft as fft
#load blur.txt
txt=np.loadtxt('/Users/Matt/Desktop/PHY407/blur.txt',unpack=True)
#transpose data to get the orientation correct
txt=np.transpose(txt)
#show image and make sure theres a gright spot as the top
plt.imshow(txt,cmap='gray')
plt.xticks([])
plt.yticks([])
plt.title('Blurry image')
plt.savefig('blurry.jpeg')
plt.show()
#b
#define a function to calculate the gaussian 
def psf(mat,sigma):    
    """This function takes a matrix and a sigma and outputs the 
    2D gaussian matrix with the given sigma onn a matrix of the 
    same size."""
    rows=len(mat)
    cols=len(mat[0])
    gauss=np.empty([rows,cols])
    for i in range(rows):
        ip=i
        if ip>rows/2:
            ip -= rows#bottom half of rows moved to negative numbers
        for j in range(cols):
            jp=j
            if jp>cols/2:
                jp -=cols#right half of columns moved to negative values
            gauss[i,j] = np.exp(-(ip**2+jp**2)/(2.*sigma**2))#compute gaussian
    return gauss
#calculate the matching gaussian function to out image
gauss=psf(txt,25)
#show corresponding gaussian matrix
plt.imshow(gauss,cmap='gray')
plt.title('Gaussian')
plt.savefig('gaussian.jpeg')
plt.show()
#c
#rfft
#use the rfft functions to fft the image and the gaussian
matfft=fft.rfft2(txt)
gssfft=fft.rfft2(gauss)
#set epsilon for the division loop
eps=1e-3
matdiv=matfft.copy()
#use a loop to divide the fft of the immahe by the fft of the gaussian
for i in range(len(gssfft)):
    for j in range(len(gssfft[0])):
        if gssfft[i][j]>eps:
            matdiv[i][j]=matfft[i][j]/gssfft[i][j]
new=fft.irfft2(matdiv)#inverse fourier transform to get back immage
plt.imshow(new,cmap='gray')
plt.xticks([])
plt.yticks([])
plt.title('Unblurred image')
plt.savefig('Unblurry.jpeg')