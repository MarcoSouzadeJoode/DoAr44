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
    frequency =  (344429.009 -64*15.625) * 1e6 * u.Hz
    
    BT = (I_nu_matrix.T).to(u.K,
                equivalencies=u.brightness_temperature(frequency)).value
    
    
    
    
    # Open the file
    with open("../ALMA_reduced.dat.vis2.syn.dat", 'r') as file:
        # Read the lines
        lines = file.readlines()[1:]
    
    # Initialize lists to store data
    ucoord = []
    vcoord = []
    hjd = []
    eff_wave = []
    vis2data = []
    vis2err = []
    vis2syn = []
    chi2 = []
    
    # Iterate over each line
    for line in lines:
        # Split the line into individual values
        values = line.split()
        
        # Convert each value to the appropriate data type and append to respective lists
        ucoord.append(float(values[0]))
        vcoord.append(float(values[1]))
        hjd.append(float(values[2]))
        eff_wave.append(float(values[3]))
        vis2data.append(float(values[4]))
        vis2err.append(float(values[5]))
        vis2syn.append(float(values[6]))
        chi2.append(float(values[7]))
    
    # Convert lists to NumPy arrays
    ucoord = np.array(ucoord)
    vcoord = np.array(vcoord)
    hjd = np.array(hjd)
    eff_wave = np.array(eff_wave)
    vis2data = np.array(vis2data)
    vis2err = np.array(vis2err)
    vis2syn = np.array(vis2syn)
    chi2 = np.array(chi2)
    
    uwave = ucoord / eff_wave
    vwave = vcoord / eff_wave
    
    CHI2 = np.sum(chi2)
    print(CHI2)
    
    uvdist = np.sqrt(uwave**2 + vwave**2) * 1e-6
    
    
 
    file_path = "../SED_NEBULA_DR_cgs_flux_inner.sed.syn.dat"

    data = {'wave': [], 'flux': [], 'fluxsyn': [], "chi2":[]}
    
    with open(file_path, "r") as f:
        lines = f.readlines()[1:]
        
    for line in lines:
        values = line.split()
        data['wave'].append(float(values[1]))
        data['flux'].append(float(values[3]))
        data['fluxsyn'].append(float(values[5]))
        data['chi2'].append(float(values[8]))
      
    chi2_sed = sum(data["chi2"])
    sed_wave = np.array(data['wave'])
    
    file_path2 = "../SED_NEBULA_DR_cgs_flux.sed.syn.dat"

    data2 = {'wave': [], 'flux': [], 'fluxsyn': [], "chi2":[]}
    
    with open(file_path2, "r") as f:
        lines = f.readlines()[1:]
        
    for line in lines:
        values = line.split()
        data2['wave'].append(float(values[1]))
        data2['flux'].append(float(values[3]))
        data2['fluxsyn'].append(float(values[5]))
        data2['chi2'].append(float(values[8]))
      
    chi2_sed2 = sum(data2["chi2"])
    sed_wave2 = np.array(data2['wave'])
    
    
    file_path3 = "../SED_NEBULA_DR_cgs_flux_mid.sed.syn.dat"

    data3 = {'wave': [], 'flux': [], 'fluxsyn': [], "chi2":[]}
    
    with open(file_path3, "r") as f:
        lines = f.readlines()[1:]
        
    for line in lines:
        values = line.split()
        data3['wave'].append(float(values[1]))
        data3['flux'].append(float(values[3]))
        data3['fluxsyn'].append(float(values[5]))
        data3['chi2'].append(float(values[8]))
      
    chi2_sed3 = sum(data3["chi2"])
    sed_wave3 = np.array(data3['wave'])
    
    
    

    
    
    fig, (ax4) = plt.subplots(nrows=1, ncols=1)
    
    
    
    fig.set_figheight(4.0)
    #fig.set_figwidth(13)
    fig.set_figwidth(4.0)

    
    minval = np.min(I_nu)
    print(minval)

    
    
    ax4.scatter(sed_wave[4:]*1e6, np.array(data["flux"])[4:]*0.1, c="k", marker="x",s=35, label="data")
    ax4.scatter(sed_wave[4:]*1e6, np.array(data["fluxsyn"])[4:]*0.1, c="xkcd:salmon", s=10, label="star+accretion")
    
    
    ax4.plot(sed_wave[4:]*1e6, np.array(data3["fluxsyn"])[4:]*0.1, c="blue", )
    ax4.scatter(sed_wave[4:]*1e6, np.array(data3["fluxsyn"])[4:]*0.1, c="blue",marker="v", label="inner disk", s=20)
    
    
    ax4.scatter(sed_wave[4:]*1e6, np.array(data2["fluxsyn"])[4:]*0.1, c="g", s=10, label="outer model", marker="s")


    ax4.plot(sed_wave[4:]*1e6, np.array(data["flux"])[4:]*0.1, c="k")
    ax4.plot(sed_wave[4:]*1e6, np.array(data["fluxsyn"])[4:]*0.1, c="xkcd:salmon")
    ax4.plot(sed_wave[4:]*1e6, np.array(data2["fluxsyn"])[4:]*0.1, c="g")
    
    #ax4.plot(sed_wave4*1e6, np.array(data4["flux"])*0.1, c="k")
    #ax4.plot(sed_wave4*1e6, np.array(data4["fluxsyn"])*0.1, c="r")
    
    
    


    TOTAL = np.array(data["fluxsyn"])[4:]*0.1 + np.array(data2["fluxsyn"])[4:]*0.1 + np.array(data3["fluxsyn"])[4:]*0.1
    #TOTAL = np.array(data["fluxsyn"])[4:]*0.1 + np.array(data3["fluxsyn"])[4:]*0.1
    ax4.plot(sed_wave[4:]*1e6, TOTAL, c="r")
    ax4.scatter(sed_wave[4:]*1e6, TOTAL, s=10, c="r",label="total model")

    
    #ax4.vlines(2010, 1e-14, 1e-4, color="k")
    ax4.vlines(350 * 1e3, 1e-14, 1e-4, color="k")
    
    
    ax4.set_ylabel("$F_\lambda$ (W $\\rm{m}^{-3}$)", fontsize=15)
    ax4.set_xlabel(f"$\lambda$ (um)", fontsize=15)
    ax4.set_xscale('log')
    ax4.set_yscale('log')
    
    ax4.tick_params(axis='both', which='major', labelsize=15)
    ax4.tick_params(axis='both', which='minor', labelsize=8)
    
    
    from matplotlib import ticker
    nticks = 18
    maj_loc = ticker.LogLocator(numticks=nticks)
    min_loc = ticker.LogLocator(subs='all', numticks=nticks)
    
    ax4.yaxis.set_major_locator(maj_loc)
    ax4.yaxis.set_minor_locator(min_loc)
    
    nticks = 18
    maj_loc = ticker.LogLocator(numticks=nticks)
    min_loc = ticker.LogLocator(subs='all', numticks=nticks)
    
    ax4.xaxis.set_major_locator(maj_loc)
    ax4.xaxis.set_minor_locator(min_loc)

    ax4.set_ylim(1e-15, 1e-5)
    ax4.set_xlim(0.3, 5000)
    ax4.legend(fontsize=12, frameon=True, loc='lower left', framealpha=1)
    #ax4.set_title(f"Spectral energy distribution")
    fig.tight_layout()
    

    import time
    #cb = fig.colorbar(mappable, fraction=0.046, pad=0.04)
    plt.savefig(f"../plot_logs/total_SED{time.time()}.png", dpi=500)
    plt.show()

if __name__ == "__main__":
    try:
        main(sys.argv[1])
    except:
        main(" ")











