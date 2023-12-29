#! /usr/bin/env python

import sys
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d

"""
Calculates the gauge cell size for different Ntotal, so that each mesocanonical simulation roughly contains around
70 particles for good statistics.
"""

T = 87.3 # K
kb = 1.38064852e-23 #m3kgs-2K-1
Nuc = 5 # 

def Column(filename, *args):
    df = pd.read_csv(filename, delim_whitespace=True)
    return [df.iloc[:, a] for a in args]

# For low pressures
fugacity, Npore_per_uc  = Column('GCMC_isotherm_reference.dat', 1, 4) # kPa
fugacity *= 1000
Npore = Npore_per_uc * Nuc # Total number of molecules within a simulation box

# interpolate to calculate the fugacity at a given number of particles within a pore
N_to_fugacity = interp1d(Npore, fugacity)

Ngauge = 70
N1 = np.arange(1, 502, 25)
N2 = np.arange(510, 750, 10)
Npore = np.concatenate([N1, N2])

# List of Ntotals
Ntotal_list = Npore + Ngauge

gauge_cell_size_list = []
for Nt in Ntotal_list:
    Np = Nt - Ngauge
    fug = N_to_fugacity(Np)
    # gauge cell box length calculation based on ideal gas
    Box_gauge = (Ngauge*kb*T/fug) ** (1./3)  # meters
    gauge_cell_size_list.append(Box_gauge*1e10)  # Angstroms



# save to a file
df = pd.DataFrame({'Ntotal': Ntotal_list, 'gauge_cell_size(A)': gauge_cell_size_list})
df.to_csv('Ntotal_vs_gauge_cell_size.dat', sep='\t', index=False, float_format='%15.3f')
