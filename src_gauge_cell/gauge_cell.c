#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdio.h>
#include <math.h>
#include "molecule.h"
#include "sample.h"
#include "thermo_baro_stats.h"
#include "rigid.h"


REAL GaugeCellSize;
REAL GaugeCellVolume;
REAL GaugeCellNumberOfParticles;
int GaugeCellParticleHistogram [10000];

FILE *FileHistogram;

/*********************************************************************************************************
 * Name       | PrintGaugeCellStatistics                                                                 *
 * ----------------------------------------------------------------------------------------------------- *
 * Function   | Prints gauge cell histogram to a file- HistogramGaugeCell.dat in the simulation directory*
 * Parameters | -                                                                                        *
 * Note       | -                                                                                        *
 * Refs.      | A.V. Neimark and A. Vishnyakov, - A Simulation Method for the Calculation of Chemical    *
 *            | Potentials in Small, Inhomogeneous, and Dense Systems , Journal of  Chemical Physics,    *
 *            | 2005,  V. 122, 234108.                                                                   *
 *********************************************************************************************************/

void PrintGaugeCellStatistics(void){

    if(GaugeCellSize == 0.0) return;

    int CurrentComponent = 0;
    int CurrentSystem =0;
    double Temp = (double)therm_baro_stats.ExternalTemperature[CurrentSystem];
    int Ntotal = Components[CurrentComponent].TotalNumberOfAdsorbateMolecules[CurrentSystem];
    double ChemicalPotentialAvg;
    int Ngauge, Nsystem;
    double ChemicalPotentialCanonical;
    REAL delta=NumberOfMoleculesHistogramSize[CurrentSystem]/NumberOfMoleculesRange[CurrentSystem];
    int r;

    //Nsystem = Components[CurrentComponent].NumberOfMolecules[CurrentSystem];
    GaugeCellVolume = CUBE(GaugeCellSize);
    // Calculate gauge cell histogram
    for(int i=0; i < 10000; i++) {
        r = i/delta;
        if (i %2 == 0 && (Ntotal - r) >= 0)
        {
            // number of particles in gauge cell = ntotal - number of particles in pore cell
            GaugeCellParticleHistogram[Ntotal - r] = NumberOfMoleculesHistogram[CurrentSystem][CurrentComponent][i];
        }
    }

    char buffer[256];
    sprintf(buffer, "HistogramGaugeCell.dat");
    FileHistogram = fopen(buffer, "w");
    fprintf(FileHistogram, "%20s %20s %20s %20s %20s \n", "Ntotal", "Npore", "Ngauge", "Ngauge_frequency", "Chempot[K]");

    for (int Nsim = 0; Nsim < Ntotal; Nsim++) {
        Ngauge = Ntotal - Nsim;
        if (GaugeCellParticleHistogram[Ngauge] !=0 && GaugeCellParticleHistogram[Ngauge-1] !=0 && Ngauge != 0) {
          // Number of particles in gauge cell is non-zero

            ChemicalPotentialCanonical = -Temp * log(GaugeCellVolume / Ngauge) +
                                             Temp * log((double) GaugeCellParticleHistogram[Ngauge] / GaugeCellParticleHistogram[Ngauge -1]);
            fprintf(FileHistogram, "%20d %20d %20d %20d %20.3lf\n", Ntotal, Nsim, Ngauge, GaugeCellParticleHistogram[Ngauge], ChemicalPotentialCanonical);
        }
    }
    fclose(FileHistogram);
}
