#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 29 21:38:14 2023

@author: marco
"""
import numpy as np
import matplotlib.pyplot as plt
import os
from astropy import constants as const
from astropy import units as u
import matplotlib as mpl
import sys

def img_reader(f):
    with open(f, "r+") as img:
        lines = img.readlines()
        N=len(lines)
        
        i = 0
        for line in lines:
            try:
                x, y, val = line.split()
                #print(x, y, val)
            except ValueError:
                i+= 1
        a = int(np.sqrt(N-i))
        
        img = np.zeros((a+1,a+1))
        XX = np.zeros((a+1,a+1))
        YY = np.zeros((a+1,a+1))
        j = 0
        k = 0
        for line in lines: 
            try:   
                x, y, z = line.split()
                img[j, k] = float(z)
                XX[j, k] = float(x)
                YY[j,k] = float(y)
    
            except ValueError:
                k = 0
                j+= 1
            k+= 1
    
           
    minX = min(XX.flatten()) * 6.68458712e-14
    maxX = max(XX.flatten()) * 6.68458712e-14
    
    return img, minX, maxX

def main(param):
    fn = '../tempif00/2Dimage_001'
    I_nu, minX, maxX = img_reader(fn)
    
    I_nu_matrix = I_nu * u.erg / (u.cm**2 * u.s * u.Hz * u.sr)
    # frequency = 1.499e14 * u.Hz
    frequency =  3.527e13 * u.Hz
    
    BT = (I_nu_matrix.T).to(u.K,
                equivalencies=u.brightness_temperature(frequency)).value
    
    
    
    
    fig, (ax1) = plt.subplots(nrows=1, ncols=1)
    
    
    
    fig.set_figheight(3)
    fig.set_figwidth(4)

    
    minval = np.min(I_nu)
    print(minval)
    
    
    cmap = mpl.colormaps.get_cmap('inferno')  # viridis is the default colormap for imshow
    cmap.set_bad(color='black')
    
    mappable = ax1.imshow(BT, cmap=cmap,
                          extent=[minX, maxX, minX, maxX],
                          vmax=5000,
                          aspect="equal")
    

    ax1.set_xlabel("X (au)", fontsize=15)
    ax1.set_ylabel("Y (au)", fontsize=15)
    #ax1.set_xticks([-60, -30, 0, 30, 60])
    #ax1.set_yticks([-60, -30, 0, 30, 60])
    ax1.set_title(f"Synthetic image of innermost disk\n Brightness Temperature (K)", pad=15)
    cbar = fig.colorbar(mappable, ax=ax1, fraction=0.046, pad=0.04)
    cbar.ax.set_ylabel('Brightness Temperature (K)', rotation=270, labelpad=20)


    
    fig.tight_layout()
    
  
    import time
    #cb = fig.colorbar(mappable, fraction=0.046, pad=0.04)
    plt.savefig(f"../plot_logs/inner{time.time()}.png", dpi=500)
    plt.show()

if __name__ == "__main__":
    try:
        main(sys.argv[1])
    except:
        main(" ")











