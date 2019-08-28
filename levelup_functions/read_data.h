#ifndef R_H
#define R_H

//structures
typedef struct
{
  double onShellClosure;
  double onSubshellClosure;
  double radioactive;
  double lowRelativeLevelDensity, highRelativeLevelDensity;

  int maxDistFromStability;
  int excludeStable;
}nuclide_rank_par; //parameters for ranking nulides

#endif
