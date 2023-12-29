#! /usr/bin/env python

import pandas as pd
from scipy import interpolate
import numpy as np
import matplotlib.pyplot as plt
fs = 15
plt.rcParams['font.size'] = '%s' % fs
pd.set_option('display.max_rows', None)

kb = 1.38064852e-23  # m3kgs-2K-1
temperature = 87.3  # K
print(f'Temperature is {temperature} K')


def Column(filename, *args):
    df = pd.read_csv(filename, delim_whitespace=True)
    return [df.iloc[:, a] for a in args]


def calculate_helmholtz_free_energy(Npore, chempot):
    """
    Helmholtz free energy F(N, V, T)
    F(N, V, T) = \sum _{i=0} ^ {N-1} chempot[i]
    Free energy at N is sum of chemical potentials from N=0 to N=N-1

    """
    # interpolate between x = N and y = mu
    f = interpolate.interp1d(Npore, chempot, fill_value="extrapolate")
    Nmax = int(np.floor(np.max(Npore)))

    sum = 0
    free_energy = []  # free energy at (N = 0 to Nmax with difference of 1 particle)
    chempot_list = []  # chempot at (N = 0 to Nmax with difference of 1 particle)

    N = np.arange(0, Nmax + 1, 1)  # number of particles from 0 to Nmax
    for i in N:
        chempot_i = f(i)
        sum += chempot_i
        chempot_list.append(chempot_i)
        free_energy.append(sum)

    df = pd.DataFrame({'N': N, 'free_energy': free_energy, 'chempot': chempot_list}, dtype=np.float128)

    return df


def calculate_grand_canonical_partition_function(Npore, chempot):
    """
    Grand Canonical Partition function \Xi(mu[N], V, T) = \Xi[N]
    Returns:

    """
    df = calculate_helmholtz_free_energy(Npore, chempot)

    GCPF = []  # grand canonical partition function
    # loop over chempot
    for Ni in df['N']:
        sum = 0
        # loop over N
        for Nj in df['N']:
            sum += np.exp((df['chempot'][Ni] * Nj - df['free_energy'][Nj]) / temperature, dtype=np.float128)
        GCPF.append(sum)

    return GCPF, df['N']


def calculate_true_grand_canonical_isotherm(Npore, chempot):

    df_free = calculate_helmholtz_free_energy(Npore, chempot)
    GCPF, _ = calculate_grand_canonical_partition_function(Npore, chempot)

    N_gce = []
    # loop over chempot
    for Ni in df_free['N']:
        sum = 0
        # loop over N
        for Nj in df_free['N']:
            sum += Nj * np.exp((df_free['chempot'][Ni] * Nj - df_free['free_energy'][Nj]) / temperature, dtype=np.float128)/GCPF[int(Ni)]
        N_gce.append(sum)

    fugacity = 1e30 * kb * temperature * np.exp(np.array(df_free['chempot'])/temperature) / 1000  # kPa
    df_out = pd.DataFrame({'#Npore': N_gce, 'fugacity[kPa]': fugacity, 'chempot[K]': df_free['chempot']})

    return df_out



if __name__ == '__main__':
    Ntotal, Npore, fugacity, chempot = Column('Mesocanonical.dat', 0, 1, 2, 3)

    # Calculate true GCE isotherm
    df_out = calculate_true_grand_canonical_isotherm(Npore, chempot)
    df_out.to_csv('True_GrandCanonical.dat',  sep='\t', index=False, float_format='%14.5g')
    print('True Grand canonical isotherm calculated')






