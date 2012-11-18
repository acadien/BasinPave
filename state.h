#ifndef STATE_H
#define STATE_H

#include "stdio.h"

//the state of the current collection of atoms
typedef struct state{
  float N; //how many points
  float **xi;//working space for optimizer
  float *x; //positions
  float *com;//center of mass
  float E; //potential energy
  int iters;//iterations for basin optimization
  float msd; //The msd (neighbor)
  float msdIdeal; //The msd with the ideal
}state;

void allocState(state*,int);
void freeState(state*);
void initState(state*,int);
void printState(state*,FILE*);
void printStateBounds(state*, FILE*);
void printStateEnergy(state*, FILE*);
void printStateVolume(state*, FILE*);
void copyState(state*,state*);
#endif
