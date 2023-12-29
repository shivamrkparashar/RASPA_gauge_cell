#! /usr/bin/env python

import pandas as pd
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
fs = 15
plt.rcParams['font.size'] = '%s' % fs

kb = 1.38064852e-23  # m3kgs-2K-1
temperature = 87.3  # K
print(f'Temperature is {temperature} K')


def Column(filename, *args):
    df = pd.read_csv(filename, delim_whitespace=True)
    return [df.iloc[:, a] for a in args]


class Linear:
    def __init__(self):
        self.fig, self.ax = plt.subplots(figsize=(6.4, 4.8), constrained_layout=True)

    def plot(self, x, y, sym, legendname, color='blue'):
        self.ax.plot(x, y, sym, label=legendname, color=color)

    def end(self, xlabel, ylabel, figname=' '):
        self.ax.set_xlabel(xlabel)
        self.ax.set_ylabel(ylabel)
        self.ax.grid()
        self.ax.legend(loc='best')
        self.ax.set_xlim(-8, 50)
        #self.ax.set_xlim(12, 35)
        #self.ax.set_ylim(-1600, -800)
        # self.ax.legend(bbox_to_anchor=(1.01, 0.5), loc='center left')
        self.fig.savefig('fig_%s.png' % figname, format='png', dpi=300, bbox_inches='tight')
        plt.show()


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

def calculate_energy_barriers(Npore, chempot, fugacity_kpa):
    """
    Deviation potential: W(N, V, T, \mu)
    The probability to sample (N, V, T) state in grand canonical ensemble at a given
    chemical potential \mu is given by deviation potential:
    W (N, V, T, \mu) = F (N, V, T) - \mu*N

    Args:
        Npore: list of Ns
        chempot: list of corresponding chemical potential
        fugacity_kpa: fugacity at which the energy surface is calculated

    Returns:

    """

    df_free = calculate_helmholtz_free_energy(Npore, chempot)
    chempot_i = temperature * np.log(fugacity_kpa*1000/(temperature*kb*1e30))  # K
    # fugacity_i = 1e30 * kb * temperature * np.exp(chempot_i/temperature)/1000  # kPa

    deviation_potential = (df_free['free_energy'] - chempot_i * df_free['N'])/temperature
    diff_deviation_potential = deviation_potential - np.min(deviation_potential)

    df_out = pd.DataFrame({'#N': df_free['N'], 'diff_deviation_potential': diff_deviation_potential})
    df_out.to_csv(f'deviation_potential_{fugacity_kpa}_kpa.dat', sep='\t', index=False, float_format='%15.6g')

    return df_free['N'], diff_deviation_potential


def calculate_canonical_work_function(Npore, chempot):
    """
    Canonical work function: WCE (N, V, T) = F (N, V, T) - \mu (N) * N
    WCE (N2) - WCE(N1): the difference in canonical work function determines the work of formation
    of the state with N2 particles from the state with N1 particles
    Args:
        Npore:
        chempot:

    Returns:
    """
    df_free = calculate_helmholtz_free_energy(Npore, chempot)
    canonical_work_function = (df_free['free_energy'] - df_free['chempot'] * df_free['N'])/temperature
    fugacity = 1e30 * kb * temperature * np.exp(df_free['chempot']/temperature)/1000  # kPa

    df_out = pd.DataFrame({'#fugacity[kPa]': fugacity, 'canonical_work_function': canonical_work_function})
    df_out.to_csv('canonical_work_function.dat', sep='\t', index=False, float_format='%15.6g')

    return fugacity, canonical_work_function


if __name__ == '__main__':
    Ntotal, Npore, fugacity, chempot = Column('Mesocanonical.dat', 0, 1, 2, 3)

    # Plot deviation_potential
    fugacity_kpa = 24.95  # kPa
    N, dE = calculate_energy_barriers(Npore, chempot, fugacity_kpa)
    A = Linear()
    A.plot(dE, N, '-', f'{fugacity_kpa} kPa', 'red')
    A.end('$W/k_BT$', 'N', 'deviation_potential')

    # Plot canonical work function
    fugacity_arr, canonical_work_function = calculate_canonical_work_function(Npore, chempot)
    B = Linear()
    B.plot(fugacity_arr, canonical_work_function, '-', 'Canonical work function', 'b')
    B.end('fugacity[kPa]', '$W_{CE}/k_BT$', 'canonical_work_function')
