#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <utils.h>

extern REAL GaugeCellSize;
extern REAL GaugeCellVolume;
extern REAL GaugeCellNumberOfParticles;
extern int GaugeCellParticleHistogram[10000];
extern FILE *FileHistogram;
void PrintGaugeCellStatistics(void);
