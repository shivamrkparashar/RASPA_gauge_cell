#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import glob

fs = 15
plt.rcParams['font.size'] = '%s' % fs


def Column(filename, *args):
    """
    Reads a column formatted file separated by spaces
    """
    df = pd.read_csv(filename, delim_whitespace=True)
    print(f'read {filename}')
    return [df.iloc[:, a] for a in args]


class LambdaPlot:
    def __init__(self):
        scale = 1.0
        figsize = (10*scale, 9*scale)
        fig, axs = plt.subplots(figsize=figsize, nrows=2, ncols=2, constrained_layout=True)

        self.ax1 = axs[0, 0]
        self.ax2 = axs[0, 1]
        self.ax3 = axs[1, 0]
        self.ax4 = axs[1, 1]
        self.fig = fig

    def plot(self, ax, x, y, sym='-o', legendname='', color='b', markerfacecolor='b'):
        ax.plot(x, y, sym, color=color, markerfacecolor=markerfacecolor, label=legendname)

    def end(self, xlabel, xlabel2, ylabel1, ylabel2, figname=' '):
        self.ax1.set_ylabel(ylabel1)
        self.ax2.set_xlabel(xlabel2)
        self.ax3.set_xlabel(xlabel)
        self.ax3.set_ylabel(ylabel2)

        self.ax1.set_title('Adsorption isotherms')
        self.ax2.set_title('Deviation potential')
        #self.ax3.set_title('Canonical work function')

        # set xlim and ylims
        self.ax1.set_xlim(12, 35)
        self.ax3.set_xlim(12, 35)
        self.ax1.set_ylim(300, 880)
        self.ax2.set_ylim(300, 880)
        self.ax2.set_xlim(-8, 50)
        #self.ax2.set_xlim(-8, 120)
        self.ax3.set_ylim(-1600, -800)

        self.ax1.legend(loc='best')
        self.ax2.legend(loc='best', ncol=2)
        self.ax3.legend(loc='best')

        # turn off labels and axis
        self.ax4.set_axis_off()
        self.ax1.xaxis.set_ticklabels([])
        self.ax2.yaxis.set_ticklabels([])

        self.ax1.grid()
        self.ax2.grid()
        self.ax3.grid()

        plt.savefig('fig_vertical_%s.png' %figname,format='png', dpi=300, bbox_inches='tight')
        plt.show()


A = LambdaPlot()

# Adsorption isotherm
# Mesocanonical
Ntotal, Npore, fugacity = Column('Mesocanonical.dat', 0, 1, 2)
A.plot(A.ax1, fugacity, Npore, '-o', 'Mesocanonical', 'g', 'g')
# True grand canonical
Npore, fugacity, chempot = Column('True_GrandCanonical.dat', 0, 1, 2)
A.plot(A.ax1, fugacity, Npore, '-', 'True grand canonical', 'b', 'b')

# Grand canonical
fugacity, Npore = Column('GCMC_Ads_isotherm.dat', 0, 1)
fugacity /= 1000
Npore *= 5
A.plot(A.ax1, fugacity, Npore, '-o', 'Grand canonical', 'k', 'orange')
fugacity, Npore = Column('GCMC_Des_isotherm.dat', 0, 1)
fugacity /= 1000
Npore *= 5
A.plot(A.ax1, fugacity, Npore, '-o', '', 'orange', 'white')


# Energy barriers
N, deviation_potential = Column('deviation_potential_24.95_kpa.dat', 0, 1)
A.plot(A.ax2, deviation_potential, N, '-', '25 kPa', 'b', 'b')

# Canonical work function
fugacity, canonical_work_function = Column('canonical_work_function.dat', 0, 1)
A.plot(A.ax3, fugacity, canonical_work_function, '-', 'Canonical work function', 'm', 'm')


A.end('Fugacity (kPa)', '$W/k_BT$',  'N (molec/sc)', '$W_{CE}/k_BT$', 'Isotherm_DeviationPotential_WorkFunction')
