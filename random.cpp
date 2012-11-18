#include "random.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "unistd.h"
#include "time.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;
using namespace Random;

void initrng(){
  //gsl_rng_env_setup();
  long seed=time(NULL)*getpid();
  T=gsl_rng_mt19937;
  randstor=gsl_rng_alloc(T);
  gsl_rng_set(randstor,seed);
}

float mrand(){//Mersenne Twister random number
  return gsl_rng_uniform(randstor);
}  

float mnormrand(float sig){//Guassian random number
  return gsl_ran_gaussian(randstor,sig);
}
