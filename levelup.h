#ifndef LU_H
#define LU_H

#include <signal.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <readline/readline.h>
#include <readline/history.h>
#include <stdint.h> //allows uint8_t and similiar types

//spectrum data file specs
#define S32K 32768
#define NSPECT 100

//increasing these numbers will increase the size of 
//the nuclear database stored in memory
#define MAXCASCDELENGTH 20
#define MAXCASCDESPERNUCL 50
#define MAXGAMMASPERLEVEL 10
#define MAXSPPERLEVEL 3
#define MAXNUMNUCL 3500
#define MAXNUMLVLS 200000

//structures
typedef struct
{
  short nucl[MAXNUMNUCL]; //nuclide indices
  float rankVal[MAXNUMNUCL]; //nuclide rank values
}nuclide_rank_info; //an individual gamma cascade

typedef struct
{
  short numLevels; //number of steps in the cascade
  float energies[MAXCASCDELENGTH]; //energies of the levels in the cascade in keV
  float gammaEnergies[MAXCASCDELENGTH];
}gamma_cascade; //an individual gamma cascade

typedef struct
{
  short halfInt; //if spin-parity value is half integer (1=true, in this case spinVal is multiplied by 0.5)
  short spinVal; //the spin value (-1 if unknown)
  short parVal;  //the parity value (1 if positive, -1 if negative, 0 if unknown)
  short tentative; //0 if not tentative, 1 if all tentative, 2 if only spin tentative, 3 if GE
}spinparval; //spin parity value

typedef struct
{
  float energy; //level energy in keV
  short energyerr; //energy uncertainty value
  float lifetime; //level lifetime (-1 if unknown)
  short lifetimeUnit; //units for level lifetime (-1=stable,0=years,1=days,2=hours,3=minutes,4=s,5=ms,6=us,7=ns,8=ps,9=fs,10=as,11=eV,12=keV,13=MeV)
  short lifetimeErr; //lifetime uncertainty value
  short numspinparvals; //number of assigned spin parity values
  spinparval spval[MAXSPPERLEVEL]; //assinged spin parity value(s) 
  short numGammas; //number of gamma rays in this level
  float gamma_energies[MAXGAMMASPERLEVEL];
  float gamma_intensities[MAXGAMMASPERLEVEL];
}level; //an individual excited level

typedef struct
{
  char nuclName[10]; //name of the nuclide, eg. '68SE'
  short N; //neutrons in nuclide
  short Z; //protons in nuclide
  float qbeta, qalpha;
  float sp, sn; //proton and neutron separation energies
  uint32_t firstLevel; //index of first level in this nuclide
  uint16_t numLevels; //number of excited levels in this nuclide
  short numCascades; //number of cascades stored for this nuclide
  gamma_cascade cascades[MAXCASCDESPERNUCL]; //cascades belonging to the nuclide
}nucl; //gamma data for a given nuclide

typedef struct
{
  int16_t numNucl; //number of nuclides for which data is stored (-1 if no nuclides)
  uint32_t numLvls; //number of levels across all nuclides
  nucl nuclData[MAXNUMNUCL]; //data for individual nuclides
  level levels[MAXNUMLVLS]; //levels belonging to nuclides
}ndata; //complete set of gamma data for all nuclides


typedef struct
{
  short numSpectra; //number of spectra
  double hist[NSPECT][S32K]; //spectrum histogram data
  double sumHist[S32K];
}inp_sp;

#endif

