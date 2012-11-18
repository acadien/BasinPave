#ifndef LOCALMIN_H
#define LOCALMIN_H

#include "state.h"

namespace Potential{
  static float xi,xj,yi,yj,zi,zj;
  static float dx,dy,dz,r,ir,ir3;
}

//Args for LJ/powell function call
typedef struct ARGST{
  int N;
  int d;
  float alpham;
  state *sp,*spp;

  //Pavement tools
  int *eBins,nBin;
  float minE,maxE,delE;
}ARGST;

float LJpot(float* cs, void* args);
float LJpotPunish(float* cs, void* args);
float LJpotAtom(float* cs, void* args);

void basinPowell(state* s,float ftol, float (*func)(float [], void*),void* args);
void basinJiggle(state* s, float (*func)(float [], void*),void* args);

#endif
