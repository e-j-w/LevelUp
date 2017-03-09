#ifndef DBO_H
#define DBO_H

typedef struct
{
	int numMatchingGammas[MAXNUMNUCL][MAXCASCDESPERNUCL];
  double matchCentroid[MAXNUMNUCL][MAXCASCDESPERNUCL][MAXCASCDELENGTH];
}casc_match_par; //fit parameters

#endif
