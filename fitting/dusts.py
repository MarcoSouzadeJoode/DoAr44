#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 15:27:47 2024

@author: marco
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  1 15:14:42 2023

@author: marco
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.patches as mp

# List of files to process
files = [
    "/home/marco/Desktop/ALMA/DoAr44_Souza_Broz/fitting/pyrmg40_opac_all.txt",
    "/home/marco/Desktop/ALMA/DoAr44_Souza_Broz/fitting/pyrmg80_opac_all.txt",
    "/home/marco/Desktop/ALMA/DoAr44_Souza_Broz/fitting/perovskite_opac_all.txt",
    "/home/marco/Desktop/ALMA/DoAr44_Souza_Broz/fitting/ammonia_opac_all.txt",
    "/home/marco/Desktop/ALMA/DoAr44_Souza_Broz/fitting/forsterite_opac_all.txt",
    "/home/marco/Desktop/ALMA/DoAr44_Souza_Broz/fitting/iron_opac_all.txt",
    "/home/marco/Desktop/ALMA/DoAr44_Souza_Broz/fitting/waterice_opac_all.txt",
]

# Function to process and plot data for a single file
def process_and_plot(file_path, figure_index):
    data = []

    with open(file_path, 'r') as file:
        current_block = []
        for line in file:
            if line.strip():
                columns = line.strip().split()
                current_block.append([float(col) for col in columns])
            else:
                if current_block:
                    block_data = {
                        'Column1': [row[0] for row in current_block],
                        'Column2': [row[1] for row in current_block],
                        'Column3': [row[2] for row in current_block],
                        'Column4': [row[3] for row in current_block],
                        'Column5': [row[4] for row in current_block],
                    }
                    data.append(block_data)
                    current_block = []

    plt.figure(figure_index)
    for idx in range(0, min(len(data), 20), 2):
        lab = data[idx]["Column1"][0]
        plt.plot(data[idx]['Column2'], data[idx]['Column4'], label=lab)

    plt.yscale("log")
    plt.xscale("log")
    plt.ylabel("Opacity")
    plt.xlabel("Wavelength ($\mu$m)")

    plt.xlim(1e-1, 1.1e3)
    B7 = mp.Rectangle([800, 1e-9], 300, 1e5, color="r", fill=True, alpha=0.3)
    B8 = mp.Rectangle([600, 1e-9], 200, 1e5, color="g", fill=True, alpha=0.3)
    B9 = mp.Rectangle([400, 1e-9], 100, 1e5, color="b", fill=True, alpha=0.3)
    B10 = mp.Rectangle([300, 1e-9], 100, 1e5, color="orange", fill=True, alpha=0.3)
    SPITZER = mp.Rectangle([5, 1e-4], 25, 1e5, color="r", fill=True, alpha=0.3)
    NIRSPEC = mp.Rectangle([0.6, 1e-4], 5-0.6, 1e5, color="b", fill=True, alpha=0.3)
    
    
    plt.gca().add_patch(B7)
    plt.gca().add_patch(B8)
    plt.gca().add_patch(B9)
    plt.gca().add_patch(B10)
    plt.gca().add_patch(SPITZER)
    plt.gca().add_patch(NIRSPEC)

    plt.vlines(870, 1e-7, 1e5, color="r")
    plt.ylim(1e-7, 1e5)
    plt.legend(ncol=2)
    plt.title(f"Opacity Plot for {file_path.split('/')[-1]}")
    plt.show()

# Loop through the files and generate plots
for i, file_path in enumerate(files):
    process_and_plot(file_path, figure_index=i+1)
