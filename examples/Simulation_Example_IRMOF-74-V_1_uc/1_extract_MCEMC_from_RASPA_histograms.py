#! /usr/bin/env python

import os
import numpy as np
import pandas as pd
from pathlib import Path
import glob

T = 87.3  # K
kb = 1.38064852e-23  # m3kgs-2K-1
Adsorbate = 'argon'


def Column(filename, *args):
    df = pd.read_csv(filename, delim_whitespace=True)
    return [df.iloc[:, a] for a in args]


def calculate_fugacity(ngauge, freq, gauge_cell_size):
    """
    Args:
        ngauge: numpy array of Number of particles in gauge cell
        freq: numpy array of corresponding frequency
        gauge_cell_size: length of the gauge cell

    Returns: Ngauge_avg, fugacity in kPa

    """
    avg_Ng = np.dot(ngauge, freq) / np.sum(freq)
    avg_chempot = -T * np.log(gauge_cell_size ** 3 / (avg_Ng))  # K
    fugacity = 1e30 * kb * T * np.exp(avg_chempot / T) / 1000  # kPa

    return avg_Ng, avg_chempot, fugacity


def MCEMC_from_RASPA_histogram(Nt, Lg, RASPA_histogram_file_path):
    """
    Calculates the mean density gauge cell (MDGC) points for a given RASPA histogram file
    Args:
        Nt: Total number of particles (Nt = Npore + Ngauge)
        Lg: Gauge cell box length in Angstroms
        RASPA_histogram_file_path: path to RASPA histogram file

    Returns: one row of dataframe with Ntotal, Npore_avg, and fugacity(kPa)
    """
    if not Path(RASPA_histogram_file_path).exists():
        print('%s not found. This simulation has not finished, skipping ' % (RASPA_histogram_file_path))
        return

    npore, freq = Column(RASPA_histogram_file_path, 0, 1)
    ngauge_avg, chempot, fugacity = calculate_fugacity(Nt - npore, freq, Lg)

    npore_avg = Nt - ngauge_avg

    # one row of dataframe
    df = pd.DataFrame({'#Ntotal': Nt, 'Npore': npore_avg, 'fugacity[kPa]': fugacity, 'chempot[K]': chempot}, index=[0])

    return df


if __name__ == '__main__':

    df = pd.DataFrame()
    Ntotal, Lgauge = Column('Ntotal_vs_gauge_cell_size.dat', 0, 1)

    for Nt, Lg in zip(Ntotal, Lgauge):
        try:
            raspa_histogram_file = glob.glob('N_%d/NumberOfMoleculesHistograms/System_0/*' % Nt)[0]
        except:
            print('Simulation Not completed, skipping'); continue

        df = pd.concat([df, MCEMC_from_RASPA_histogram(Nt, Lg, raspa_histogram_file)])

    df = df.sort_values(by='#Ntotal', ignore_index=True)
    df.to_csv('Mesocanonical.dat', sep='\t', index=False, float_format='%1.3f')
    print('Mesocanonical.dat printed')
