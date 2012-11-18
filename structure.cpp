#include "math.h"
#include "structure.h"
#include "state.h"
#include "constants.h"
#include "localMin.h"
#include "random.h"

using namespace std;
using namespace Struc;

//Distance between 2 points
float dist(float *a, float* b){
  strucx=a[0]-b[0];
  strucy=a[1]-b[1];
  strucz=a[2]-b[2];
  return sqrt(strucx*strucx + strucy*strucy + strucz*strucz);
}

//Distance between point and origin
float origDist(float *a){
  return sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
}

//Calculate the com
void com(state* s){
  //First re-center about the center of mass.
  float com[3]={0.0,0.0,0.0};
  for(int i=0;i<s->N;i++){
    com[0]+=s->x[3*i+1];
    com[1]+=s->x[3*i+2];
    com[2]+=s->x[3*i+3];
  }
  for(int i=0;i<3;i++)
    s->com[i]=com[i]/s->N;
}

void recenter(state* s){
  com(s);
  float *com=s->com;
  for(int i=0;i<s->N;i++){
    s->x[3*i+1]-=com[0];
    s->x[3*i+2]-=com[1];
    s->x[3*i+3]-=com[2];
  }
  for(int i=0;i<3;i++)
    s->com[i]=0.;
}

float msd(state* a, state* b){
  float d=0,m=0;
  float* xa=a->x;
  float* xb=b->x;
  int N=a->N;
  for(int i=0;i<N;i++){
    d=dist(&xa[3*i+1],&xb[3*i+1]);
    m+=d*d;
  }
  return sqrt(m/(float)N);
}

//Move atoms from outside the sphere... inside sphere somewhere
float salt(state* s){
  int N = s->N,cnt,outcount=0;
  float *x=s->x;
  float r,theta,phi;
  ARGST args;
  args.N=N;

  for(int i=0;i<N;i++){
    if( origDist(&(x[3*i+1])) > boundr ){
      outcount++;
      cnt=0;
      args.d=i;
      while(true){
	cnt++;
	r=boundr*mnormrand(1.0);
	while(r>boundr) r-=boundr;
	while(r<-boundr) r+=boundr;
	theta=2*M_PI*mrand();
	phi=acos(2.0*mrand()-1.0);
	x[i*3+1]=r*sin(theta)*cos(phi);
	x[i*3+2]=r*sin(theta)*sin(phi);
	x[i*3+3]=r*cos(theta);
	int a=LJpot(x,(void*)&args);
	if(a - s->E < 5.)
	  break;
	else
	  printf("???? salt fail %f ????\n",a-s->E);
	if(cnt>10000){
	  printf("Error:unable to salt atom!\n");
	  exit(0);
	}
      }
    }
  }
  if(outcount>0)
    printf("salt cnt=%d\n",outcount);
}

//compresses axes that are greater than the average boundary width
//intended to be used during initialization
void cubify(state *s){
  mnx=1E10;
  mny=1E10;
  mnz=1E10;
  mxx=-1E10;
  mxy=-1E10;
  mxz=-1E10;
  float *xs=s->x;
  for(int i=0;i<s->N;i++){
    if(xs[3*i+1]<mnx)
      mnx=xs[3*i+1];
    if(xs[3*i+1]>mxx)
      mxx=xs[3*i+1];
    if(xs[3*i+2]<mny)
      mny=xs[3*i+2];
    if(xs[3*i+2]>mxy)
      mxy=xs[3*i+2];
    if(xs[3*i+3]<mnz)
      mnz=xs[3*i+3];
    if(xs[3*i+3]>mxz)
      mxz=xs[3*i+3];
  }
  strucx=mxx-mnx;
  strucy=mxy-mny;
  strucz=mxz-mnz;
  float avgwidth=(strucx+strucy+strucz)/3.0;
  strucx=0.5*(1.+avgwidth/strucx);
  strucy=0.5*(1.+avgwidth/strucy);
  strucz=0.5*(1.+avgwidth/strucz);
  for(int i=0;i<s->N;i++){
    xs[3*i+1]*=strucx;
    xs[3*i+2]*=strucy;
    xs[3*i+3]*=strucz;
  }
}
