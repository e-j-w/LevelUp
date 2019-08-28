#ifndef LU_H
#define LU_H

#include <signal.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

//spectrum data file specs
#define S32K 32768
#define NSPECT 100

//increasing these numbers will increase the size of 
//the nuclear database stored in memory
#define MAXCASCDELENGTH 20
#define MAXCASCDESPERNUCL 50
#define MAXGAMMASPERLEVEL 10
#define MAXSPPERLEVEL 3
#define MAXLEVELSPERNUCL 200
#define MAXNUMNUCL 3500

//structures
typedef struct
{
  int numLevels; //number of steps in the cascade
  double energies[MAXCASCDELENGTH]; //energies of the levels in the cascade in keV
  double gammaEnergies[MAXCASCDELENGTH];
}gamma_cascade; //an individual gamma cascade

typedef struct
{
  int halfInt; //if spin-parity value is half integer (1=true, in this case spinVal is multiplied by 0.5)
  int spinVal; //the spin value (-1 if unknown)
  int parVal;  //the parity value (1 if positive, -1 if negative, 0 if unknown)
  int tentative; //0 if not tentative, 1 if all tentative, 2 if only spin tentative, 3 if GE
}spinparval; //spin parity value

typedef struct
{
  double energy; //level energy in keV
  int energyerr; //energy uncertainty value
  double lifetime; //level lifetime (-1 if unknown)
  int lifetimeUnit; //units for level lifetime (-1=stable,0=years,1=days,2=hours,3=minutes,4=s,5=ms,6=us,7=ns,8=ps,9=fs,10=as,11=eV,12=keV,13=MeV)
  int lifetimeErr; //lifetime uncertainty value
  int numspinparvals; //number of assigned spin parity values
  spinparval spval[MAXSPPERLEVEL]; //assinged spin parity value(s) 
  int numGammas; //number of gamma rays in this level
  double gamma_energies[MAXGAMMASPERLEVEL];
  double gamma_intensities[MAXGAMMASPERLEVEL];
}level; //an individual excited level

typedef struct
{
	char nuclName[10]; //name of the nucleus, eg. '68SE'
	int N; //neutrons in nucleus
	int Z; //protons in nucleus
	int numLevels; //number of excited levels in this nucleus
	level levels[MAXLEVELSPERNUCL]; //levels belonging to the nucleus
	int numCascades; //number of cascades stored for this nucleus
	gamma_cascade cascades[MAXCASCDESPERNUCL]; //cascades belonging to the nucleus
}nucl; //gamma data for a given nucleus

typedef struct
{
	int numNucl; //number of nuclei for which data is stored
	nucl nuclData[MAXNUMNUCL];
}ndata; //complete set of gamma data for all nuclei


typedef struct
{
	int numSpectra; //number of spectra
	double hist[NSPECT][S32K]; //spectrum histogram data
	double sumHist[S32K];
}inp_sp;

#endif

