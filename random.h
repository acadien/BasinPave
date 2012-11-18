#ifndef RANDOM_H
#define RANDOM_H

#include "math.h"
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

namespace Random{
  static const gsl_rng_type* T;
  static gsl_rng* randstor;
}

void initrng();
float mrand();//Mersenne Twister random number
float mnormrand(float);//Guassian random number

#endif
