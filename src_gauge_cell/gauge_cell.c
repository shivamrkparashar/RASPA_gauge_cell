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
    int Ngauge, Nsystem;
    Nsystem = Components[CurrentComponent].NumberOfMolecules[CurrentSystem];
    float ChemicalPotentialCanonical;

    REAL delta=NumberOfMoleculesHistogramSize[CurrentSystem]/NumberOfMoleculesRange[CurrentSystem];
    int r;
    for(int i=0; i < 1000; i++) {
        r = i/delta;
        //ParticleHistogram[r] = NumberOfMoleculesHistogram[CurrentSystem][CurrentComponent][i];
        //GaugeCellParticleHistogram[Ntotal - r] = ParticleHistogram[r];
        if (i %2 == 0 && (Ntotal - r) >= 0)
        {
            GaugeCellParticleHistogram[Ntotal - r] = NumberOfMoleculesHistogram[CurrentSystem][CurrentComponent][i];
            //printf("%d  %20lf  %20d \n", i,  NumberOfMoleculesHistogram[CurrentSystem][CurrentComponent][i], r) ;
            //printf("%20lf  %20d \n", (double)r, GaugeCellParticleHistogram[Ntotal - r]);
        }
        //if(NumberOfMoleculesHistogram[CurrentSystem][CurrentComponent][i]>0.0)
        //    printf("%20lf  %20lf \n", (double)r, NumberOfMoleculesHistogram[CurrentSystem][CurrentComponent][i]);
        //
        //if (i <= Ntotal*delta)
        //    printf("%d  %20lf  %20d \n", i,  NumberOfMoleculesHistogram[CurrentSystem][CurrentComponent][i], r) ;
        
    }
    // Print results on screen
    
    /*
    for(int r=0; r < Ntotal; r++){
            
            printf("%d  %20d \n", r, GaugeCellParticleHistogram[r]) ;
            
            }
            */

    char buffer[256];
    printf("%20s %20s %20s %20s %20s \n", "MaxNparticles", "Npore", "Ngauge", "Ngauge_frequencey", "Chempot");
    sprintf(buffer, "HistogramGaugeCell.dat");
    FileHistogram = fopen(buffer, "w");

    fprintf(FileHistogram, "%20s %20s %20s %20s %20s \n", "MaxNparticles", "Npore", "Ngauge", "Ngauge_frequencey", "Chempot");
    //printf("%20d %20d %20d \n", Ntotal, Nsystem, Ntotal-Nsystem);


    for (int Nsim = 0; Nsim < Ntotal; Nsim++) {
        Ngauge = Ntotal - Nsim;
        if (GaugeCellParticleHistogram[Ngauge] !=0 && GaugeCellParticleHistogram[Ngauge-1] !=0 && Ngauge != 0) {
          // Number of particles in gauge cell is non-zero

            ChemicalPotentialCanonical = -Temp * log(GaugeCellVolume / Ngauge) +
                                             Temp * log((double) GaugeCellParticleHistogram[Ngauge] / GaugeCellParticleHistogram[Ngauge -1]);
            fprintf(FileHistogram, "%20d %20d %20d %20d %20.3lf\n", Ntotal, Nsim, Ngauge, GaugeCellParticleHistogram[Ngauge], ChemicalPotentialCanonical);
            printf("%20d %20d %20d %20d %20.3lf\n", Ntotal, Nsim, Ngauge, GaugeCellParticleHistogram[Ngauge], ChemicalPotentialCanonical);
        }
    }
    fclose(FileHistogram);
}
//ChemicalPotentialAvg = -Temp * log(GaugeCellVolume / Lambda3_Sigma3 / (Ntotal-AverageNumberOfParticle + 1) );
