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
    fn = '../tempspe00/2Dimage_001'
    fn="../tempif00/2Dimage_001"
    I_nu, minX, maxX = img_reader(fn)
    
    I_nu_matrix = I_nu * u.erg / (u.cm**2 * u.s * u.Hz * u.sr)
    #frequency = 1.499e14 * u.Hz
    frequency = (344429.009 -64*15.625) * 1e6 * u.Hz
    
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
    #file_path = "SED_HOT.sed.syn.dat"
    #file_path = "/home/marco/Desktop/ALMA/oct2023/Nebula2/fitting/SED_STAR_ONLY_DR.dat.sed.syn.dat"
    
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
    
    
    
    
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2)
    
    
    
    fig.set_figheight(4)
    fig.set_figwidth(8)
    
    
    minval = np.min(I_nu)
    print(minval)
    
    
    cmap = mpl.colormaps.get_cmap('viridis')  # viridis is the default colormap for imshow
    cmap.set_bad(color='black')
    
    mappable = ax1.imshow((BT), cmap=cmap,
                          extent=[minX, maxX, minX, maxX],
                          vmax = 50,
                          aspect="equal")
    
    
    ax1.set_xlabel("X (au)", fontsize=10)
    ax1.set_ylabel("Y (au)", fontsize=10)
    #ax1.set_xticks([-60, -30, 0, 30, 60])
    #ax1.set_yticks([-60, -30, 0, 30, 60])
    ax1.set_title(f"Synthetic shellspec image \n Brightness Temperature (K)")
    fig.colorbar(mappable, ax=ax1,fraction=0.046, pad=0.04)
    
    ax2.scatter(uvdist, vis2data, c="k", s=1,alpha=0.5, label="ALMA Band 7")
    ax2.scatter(uvdist, vis2syn, c="r", s=3, label="synthetic")
    #ax2.legend(fontsize=10)
    ax2.set_ylim(-0.1, 1.1)
    ax2.set_xlim(0, 2)
    
    ax2.set_ylabel("Visibility $V^2$", fontsize=10)
    ax2.set_xlabel("UV distance $(M\lambda)$", fontsize=10)
    ax2.set_title(f"Interferometric visibilities \n $\chi^2 = $ {CHI2:.3e}")
    ax2.legend()
    
    

    """
    
    ax4.scatter(sed_wave*1e9, np.array(data["flux"])*0.1, c="k", s=3, label="data")
    ax4.scatter(sed_wave*1e9, np.array(data["fluxsyn"])*0.1, c="r", s=3, label="synthetic")
    
    
    #ax4.vlines(2010, 1e-14, 1e-4, color="k")
    ax4.vlines(656, 1e-8, 1e-4, color="k")
    
    
    ax4.set_ylabel("$F_\lambda$ (W $\\rm{m}^3$)", fontsize=10)
    ax4.set_xlabel(f"$\lambda$ (nm)", fontsize=10)
    ax4.set_xscale('log')
    ax4.set_yscale('log')
    ax4.legend(fontsize=10, frameon=True)
    ax4.set_title(f"Spectral energy distribution \n $\chi^2 = $ {chi2_sed:.1e}")
    """
    fig.tight_layout()
    
    
    import time
    #cb = fig.colorbar(mappable, fraction=0.046, pad=0.04)
    plt.savefig(f"../plot_logs/SYN_IF_SPE_{time.time()}.png", dpi=500)
    plt.show()

if __name__ == "__main__":
    try:
        main(sys.argv[1])
    except:
        main(" ")










