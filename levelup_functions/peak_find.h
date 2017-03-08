#ifndef PEAK_FIND
#define PEAK_FIND

#include <stdlib.h>
#include <stdio.h>

#define BIG_NUMBER 1E10
#define MAXPEAKSTOFIND 100

typedef struct
{
	double filterValue[S32K];
	int numPeaksFound;
  int centroid[MAXPEAKSTOFIND];
}peak_fit_par; //fit parameters


peak_fit_par findPeak(const double*,double);

#endif
