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
int GaugeCellParticleHistogram [1000];

FILE *FileHistogram;

void PrintGaugeCellStatistics(void){

    if(GaugeCellSize == 0.0) return 

    //GaugeCellVolume = GaugeCellSize*GaugeCellSize*GaugeCellSize;
    GaugeCellVolume = CUBE(GaugeCellSize);
    int CurrentComponent = 0;
    int CurrentSystem =0;
    double Temp = (double)therm_baro_stats.ExternalTemperature[CurrentSystem];
          
    int Ntotal = Components[CurrentComponent].CreateNumberOfMolecules[CurrentSystem];
    double ChemicalPotentialAvg;

    for(int i=0; i < 1000; i++) {
    GaugeCellParticleHistogram[i] = NumberOfMoleculesHistogram[CurrentSystem][CurrentComponent][i];
    }
    // Print results on screen
    int Ngauge, Nsystem;
    Nsystem = Components[CurrentComponent].NumberOfMolecules[CurrentSystem];
    float ChemicalPotentialCanonical;

    char buffer[256];
    printf("%20s %20s %20s %20s \n", "MaxNparticles", "Npore", "Ngauge", "Chempot");
    sprintf(buffer, "HistogramGaugeCell.dat");
    FileHistogram = fopen(buffer, "w");

    fprintf(FileHistogram, "%20s %20s %20s %20s \n", "MaxNparticles", "Npore", "Ngauge", "Chempot");
    printf("%20d %20d %20d \n", Ntotal, Nsystem, Ntotal-Nsystem);


    for (int Nsim = 0; Nsim < Ntotal; Nsim++) {
        Ngauge = Ntotal - Nsim;
        if (GaugeCellParticleHistogram[Ngauge-1] !=0 && GaugeCellParticleHistogram[Ngauge-2] !=0) {
          // Number of particles in gauge cell is non-zero

            ChemicalPotentialCanonical = -Temp * log(GaugeCellVolume / Ngauge) +
                                             Temp * log((double) GaugeCellParticleHistogram[Ngauge-1] / GaugeCellParticleHistogram[Ngauge - 2]);
            fprintf(FileHistogram, "%20d %20d %20d %20.3lf\n", Ntotal, Nsim, Ntotal - Nsim, ChemicalPotentialCanonical);
            printf("%20d %20d %20d %20.3lf\n", Ntotal, Nsim, Ntotal - Nsim, ChemicalPotentialCanonical);
        }
    }
    fclose(FileHistogram);
}
//ChemicalPotentialAvg = -Temp * log(GaugeCellVolume / Lambda3_Sigma3 / (Ntotal-AverageNumberOfParticle + 1) );
