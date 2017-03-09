#ifndef PEAK_FIND
#define PEAK_FIND

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define BIG_NUMBER 1E10
#define MAXPEAKSTOFIND 100

typedef struct
{
	double filterValue[S32K],filter2Value[S32K];
	int numPeaksFound;
  double centroid[MAXPEAKSTOFIND];
  double intensity[MAXPEAKSTOFIND];
}peak_fit_par; //fit parameters

#endif
