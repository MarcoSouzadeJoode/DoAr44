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
    
    
    
    
    # Open the file
    with open("MATISSE_N2.vis2.syn.dat", 'r') as file:
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
    if_eff_wave = np.array(eff_wave)
    vis2data = np.array(vis2data)
    vis2err = np.array(vis2err)
    vis2syn = np.array(vis2syn)
    chi2 = np.array(chi2)
    
    uwave = ucoord / eff_wave
    vwave = vcoord / eff_wave
    
    CHI2 = np.sum(chi2)
    print(CHI2)
    
    uvdist = np.sqrt(uwave**2 + vwave**2) * 1e-6
    
    
    fn = "XSHOOTER_HALFA.spe.syn.dat"
    #fn = "/home/marco/Desktop/ALMA/oct2023/Nebula2/fitting/UVES_HBETA.spe.syn.dat"
    with open(fn, 'r') as file:
        # Read the lines
        lines = file.readlines()[1:]
    
    # Initialize lists to store data
    hjd2 = []
    eff_wave = []
    eff_band = []
    flux = []
    error = []
    fluxsyn = []
    dataset = []
    phase = []
    chi2_spe = []
    
    # Iterate over each line
    for line in lines:
        # Split the line into individual values
        values = line.split()
        
        # Convert each value to the appropriate data type and append to respective lists
        hjd2.append(float(values[0]))
        eff_wave.append(float(values[1]))
        eff_band.append(float(values[2]))
        flux.append(float(values[3]))
        error.append(float(values[4]))
        fluxsyn.append(float(values[5]))
        dataset.append(float(values[6]))
        phase.append(float(values[7]))
        chi2_spe.append(float(values[8]))
    
    
    
    chi22 = sum(chi2_spe)
    eff_wave = np.array(eff_wave)
    
    
    file_path = "SED_NEBULA_DR_cgs_flux_inner.sed.syn.dat"
    file_path = "SED_NEBULA_DR_cgs_outer.sed.syn.dat"

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
    
    file_path2 = "SED_NEBULA_DR_cgs_flux.sed.syn.dat"

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
    
    
    file_path3 = "SED_NEBULA_DR_cgs_flux_mid.sed.syn.dat"

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
    
    
    
    file_path4 = "SPITZER_water.sed.syn.dat"

    data4= {'wave': [], 'flux': [], 'fluxsyn': [], "chi2":[]}
    
    with open(file_path4, "r") as f:
        lines = f.readlines()[1:]
        
    for line in lines:
        values = line.split()
        data4['wave'].append(float(values[1]))
        data4['flux'].append(float(values[3]))
        data4['fluxsyn'].append(float(values[5]))
        data4['chi2'].append(float(values[8]))
      
    chi2_sed4 = sum(data4["chi2"])
    sed_wave4 = np.array(data4['wave'])
    
    
    fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3)
    
    
    
    fig.set_figheight(3)
    #fig.set_figwidth(13)
    fig.set_figwidth(10)

    
    minval = np.min(I_nu)
    print(minval)
    
    
    cmap = mpl.colormaps.get_cmap('inferno')  # viridis is the default colormap for imshow
    cmap.set_bad(color='black')
    
    mappable = ax1.imshow(BT, cmap=cmap,
                          extent=[minX, maxX, minX, maxX],
                          vmax = 5,
                          aspect="equal")
    
    
    ax1.set_xlabel("X (au)", fontsize=15)
    ax1.set_ylabel("Y (au)", fontsize=15)
    #ax1.set_xticks([-60, -30, 0, 30, 60])
    #ax1.set_yticks([-60, -30, 0, 30, 60])
    ax1.set_title(f"synthetic image @ 10 μm \n Brightness Temperature (K)", pad=15)
    fig.colorbar(mappable, ax=ax1, fraction=0.046, pad=0.04)
    
    #ax2.scatter(uvdist, vis2data, c="k", s=3,alpha=0.1, label="ALMA Band 7")
    ax2.scatter(uvdist, vis2syn, c="k", s=30, label="MATISSE N2")
    ax2.legend(fontsize=12)
    #ax2.hlines(0, 0, 2, color="b", zorder=0, linewidth=2)
    ax2.set_ylim(-0.2, 1.2)
    ax2.set_xlim(0, 18)
    
    ax2.set_ylabel("Visibility $V^2$", fontsize=15)
    ax2.set_xlabel("UV distance $(M\lambda)$", fontsize=15)
    ax2.set_title(f"Interferometric visibilities")
    
    ax3.scatter(if_eff_wave * 1e6, vis2syn, c="k", s=30, label="MATISSE N1")
    ax3.set_ylabel("Visibility $V^2$", fontsize=15)
    ax3.set_xlabel("$\lambda$ (μm)", fontsize=15)

    
    
    """
    factor = 1
    ax3.scatter(eff_wave*1e9, flux, c="k", s=3, label="VLT/XSHOOTER")
    ax3.plot(eff_wave*1e9, flux, c="k", linewidth=1, label="VLT/XSHOOTER")

    ax3.scatter(eff_wave*1e9, (np.array(fluxsyn) - 1) * factor + 1, c="r", s=3, label=f"synthetic")
    ax3.plot(eff_wave*1e9, (np.array(fluxsyn) - 1) * factor + 1, c="r", linewidth=1, label=f"synthetic")

    ax3.set_ylabel("Relative intensity", fontsize=10)
    ax3.set_xlabel(f"$\lambda$ (nm)", fontsize=10)
    
    
    
    #ax3.legend(fontsize=10)
    ax3.set_title(f"H$\\alpha$ profile \n $\chi^2 = $ {chi22:.1e}")
    ax3.set_ylim(-4, 15)
    """
    
    """
    ax4.scatter(sed_wave[:]*1e9, np.array(data["flux"])[:]*0.1, c="k", s=10, label="data")
    ax4.scatter(sed_wave[:]*1e9, np.array(data["fluxsyn"])[:]*0.1, c="g", s=10, label="inner model")
  
    
    #ax4.vlines(2010, 1e-14, 1e-4, color="k")
    ax4.vlines(350 * 1e3, 1e-14, 1e-4, color="k")
    
    
    ax4.set_ylabel("$F_\lambda$ (W $\\rm{m}^3$)", fontsize=10)
    ax4.set_xlabel(f"$\lambda$ (nm)", fontsize=10)
    ax4.set_xscale('log')
    ax4.set_yscale('log')
    ax4.set_ylim(1e-16, 1e-5)
    ax4.set_xlim(500, 5000000)
    ax4.legend(fontsize=10, frameon=True, loc='center left', bbox_to_anchor=(1, 0.5))
    ax4.set_title(f"Spectral energy distribution")
    """
    fig.tight_layout()
    
    """
    axB.scatter(uvdist, vis2data - vis2syn, c="k", alpha=1, s=0.1)
    
    
    OBS = np.array(data["flux"])[4:]*0.1
    SYN_tot = np.array(data2["fluxsyn"])[4:]*0.1 + np.array(data["fluxsyn"])[4:]*0.1
    DIFF = OBS - SYN_tot
    
    axC.set_xlim(1e5, 1e7)
    axC.set_ylim(-0.6, 0.3)

    axC.scatter(sed_wave[4:]*1e9, DIFF / OBS)
    axC.hlines(0, 2e2, 1e7, color="r")
    #axC.set_yscale("log")
    axC.set_xscale("log")
    """
    import time
    #cb = fig.colorbar(mappable, fraction=0.046, pad=0.04)
    plt.savefig(f"plots/outer{time.time()}.png", dpi=500)
    plt.show()

if __name__ == "__main__":
    try:
        main(sys.argv[1])
    except:
        main(" ")











