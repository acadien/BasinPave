#ifndef STRUCTURE_H
#define STRUCTURE_H

//Local variables, save on allocation time...
#include "state.h"

namespace Struc{
  static float strucx,strucy,strucz;
  static float mnx,mxx,mny,mxy,mnz,mxz;
}

float dist(float *a, float* b);
float origDist(float *a);
float msd(state* a, state* b);
void com(state* s);
void recenter(state* s);
float salt(state* s);
void cubify(state *s);

#endif
