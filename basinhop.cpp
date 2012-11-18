#include <stdlib.h>
#include <stdio.h>
#include <queue>
#include <string>
#include "gsl/gsl_math.h"
#include <string.h>
//mine
#include "constants.h"
#include "random.h"
#include "state.h"
#include "localMin.h"
#include "structure.h"
#include "rmsd.h"
#include "basinhop.h"

using namespace std;

int main(int argc, char *argv[]){
  state s;
  state sprime;
  state sideal;
  state* basins;

  //Settings
  int aLen=100;  //how big should the window average be
  int nAtom=atoi(argv[1]); //how many atoms
  bool MethodA;
  MethodA=false;
  char A[]="A";
  if(!strcmp(argv[2],A))
    MethodA=true;
  
  //Basin
  float ftol=0.01; //set the tolerance on basin finding algo methodA

  //MC
  int initLoop=500, hopLoop=5000; //MC loop lengths
  float MCT=0.8;                  //Monte Carlo Temperature
  float MCalpha=0.20;             //Monte Carlo initial jump length
  
  //Initialize random numbers (mersenne twist)
  initrng();

  //Load up ideal cluster
  char buf[5];
  sprintf(buf,"%d",nAtom);
  string idealfilename=buf;
  idealfilename.append("_ideal.dat");
  FILE *idealfp;
  idealfp = fopen((char*)idealfilename.c_str(),"r");
  initState(&sideal,nAtom);
  loadIdeal(&sideal,idealfp);

  //Initialize cluster
  initState(&s,nAtom);
  initState(&sprime,nAtom);
  basins=(state*)malloc(sizeof(state)*hopLoop);
  for(int i=0;i<hopLoop;i++)
    allocState(&(basins[i]),nAtom);
  ARGST *args=(ARGST*)malloc(sizeof(ARGST));
  args->N=s.N;
  args->alpham=0.2;
  args->sp=(state*)malloc(sizeof(state));
  args->spp=(state*)malloc(sizeof(state));
  initState(args->sp,nAtom);
  initState(args->spp,nAtom);

  //Pavement tools
  args->nBin=10000;
  args->eBins=(int*)malloc(sizeof(int)*args->nBin);
  args->maxE=-2.0*nAtom;
  args->minE=-5.65*nAtom;
  args->delE=(args->maxE-args->minE)/(float)args->nBin;

  //Initalize logging
  s.E=LJpotPunish(s.x,(void*)args);
  FILE* fp=stdout;
  FILE *logf;
  if(MethodA){
    if(nAtom==38)
      logf = fopen("finalPAVE_N38_A.dat","w");
    if(nAtom==76)
      logf = fopen("finalPAVE_N76_A.dat","w");
    if(nAtom==104)
      logf = fopen("finalPAVE_N104_A.dat","w");
  }
  else{
    if(nAtom==38)
      logf = fopen("finalPAVE_N38_B.dat","w");
    if(nAtom==76)
      logf = fopen("finalPAVE_N76_B.dat","w");
    if(nAtom==104)
      logf = fopen("finalPAVE_N104_B.dat","w");
  }


  //Initialize acceptance queue for dynamically altering variation parameter
  std::queue<int> accepts;
  float acceptAvg=0.5;
  for(int i=0;i<aLen/2;i++){
    accepts.push(1);
    accepts.push(0);
  }

  //Prepare optimizers
  printf("Initial\n");
  printStateVolume(&s,fp);
  printStateEnergy(&s,fp);
  printf("\n");

  //Initial cluster will be horrible, run it through some MC to reduce total energy
  for(int i=0;i<initLoop;i++){
    
    MCstep(&s,&sprime,(void*)args,ftol,MCT,MCalpha,&accepts,&acceptAvg,MethodA,true);

    printStateVolume(&s,fp);
    printStateEnergy(&s,fp);
    printStateBounds(&s,fp);
  }
  FILE* pavelog0;
  pavelog0 = fopen("pavelog0.dat","w");
  FILE* pavelog1;
  pavelog1 = fopen("pavelog1.dat","w");
  printPave(pavelog0,args);

  printf("Reduced Initial.\n");
  printStateVolume(&s,fp);
  printStateEnergy(&s,fp);
  printStateBounds(&s,fp);
  printf("\n");
  
  //Let the basin hopping begin
  float msdnow;
  for(int i=0;i<hopLoop;i++){

    MCstep(&s,&sprime,(void*)args,ftol,MCT,MCalpha,&accepts,&acceptAvg,MethodA,true);

    copyState(&s,&(basins[i]));

    printf("****************************\n");
    printf("%d\n",i);
    printStateVolume(&s,fp);
    printStateEnergy(&s,fp);
    printStateBounds(&s,fp);
    printf("****************************\n");
  }
  printPave(pavelog1,args);
  fflush(pavelog0);
  fflush(pavelog1);
  fclose(pavelog0);
  fclose(pavelog1);

  //Reoptimize with higher accuracy
  ftol=1e-3;
  for(int i=0;i<hopLoop;i++){
    printf("%d\n",i);

    copyState(&(basins[i]),&sprime);

    //if(MethodA)     //Method A
    basinPowell(&sprime,ftol,LJpot,(void*)args);
    //else            //Method B
    //  basinJiggle(&sprime,LJpot,(void*)args);

    if( sprime.E < basins[i].E )
      copyState(&sprime,&(basins[i]));
    if(i>0)
      basins[i].msd=msd(&basins[i],&basins[i-1]);

    basins[i].msdIdeal=rmsd(3*nAtom,&(basins[i].x[1]),&(sideal.x[1]));

    //Write relevant information to log
    printState(&(basins[i]),logf);
  }
  

  for(int i=0;i<hopLoop;i++)
    freeState(&(basins[i]));
  free(basins);
  if(MethodA)
    freeState(&sprime);
  freeState(&sideal);
  freeState(&s);
  freeState(args->sp);
  freeState(args->spp);
  free(args->sp);
  free(args->spp);
  free(args);
  return 0;
}

void MCstep(state* s, state* sprime,void* args,float ftol, float& MCT, float& MCalpha, std::queue<int>* accepts, float* acceptAvg, bool MethodA, bool silent){

  int cnt=0;
  float alphaStep=0.0002,alphaRatio=0.99;
  float E;

  ARGST* pave=(ARGST*)args;
  float minE=pave->minE;
  float delE=pave->delE;
  int eIndex,nBin=pave->nBin;
  int N=sprime->N;
  cubify(s); //Squeeze the structure to be more "cube-like"   
  while(true){
    cnt++;

    copyState(s,sprime);

    //Salt atoms that are outside of sphere boundary
    salt(sprime);

    //Step out of local minimum!
    for(int i=0;i<N;i++){
      sprime->x[3*i+1] += (mrand()-0.5)*2.0*MCalpha;
      sprime->x[3*i+2] += (mrand()-0.5)*2.0*MCalpha;
      sprime->x[3*i+3] += (mrand()-0.5)*2.0*MCalpha;
    }
    if(MethodA)
      basinPowell(sprime,ftol,LJpot,args);
    else
      basinJiggle(sprime,LJpot,args);
    sprime->iters=0;

    //Evaluate paving
    eIndex=(int)((sprime->E-minE)/delE);
    if(eIndex>=0 and eIndex<nBin){
      sprime->E+=(float)pave->eBins[eIndex]*delE;
      printf("paveval=%d\n",pave->eBins[eIndex]);
    }
    else
      printf("----- ERROR PAVING index=%d -----\n",eIndex);

    //Calculate the Metropolis Criterion
    float weight=exp( -(sprime->E - s->E) / MCT );
    //if(!silent)
    //printf("old:%4.4f new:%4.4f | expdelE=%4.4f\n",s->E,sprime->E,weight);

    //Monte-Carlo action bam-pow
    bool accept=false;
    if(sprime->E < s->E){
      copyState(sprime,s);
      accept=true;
    }else 
      if(mrand()<weight){
	copyState(sprime,s);
	accept=true;
      }else{
	if(!silent)
	  printf("higher energy didn't take\n");
      }

    if(!MethodA and cnt>10 and sprime->E - s->E < 50){
      copyState(sprime,s);
      //basinPowell(sprime,1.0,LJpot,args);
      accept=true;
      printf("bleh\n");
    }

    //Update the acceptance queue and average acceptance
    *acceptAvg-=accepts->front()/(float)accepts->size();
    accepts->pop();
    if(accept){
      accepts->push(1);
      *acceptAvg+=1.0/(float)accepts->size();
    }
    else
      accepts->push(0);

    //Update the Monte-Carlo Temperature according to acceptance
    if(accept)
      MCalpha/=alphaRatio;    
    else if (MCalpha > alphaStep)
      MCalpha*=alphaRatio;

    if(accept){
      if(eIndex>2){
	pave->eBins[eIndex-2]+=1;
	pave->eBins[eIndex-1]+=3;
      }
      pave->eBins[eIndex]+=10;
      if(eIndex<nBin-2){
	pave->eBins[eIndex+1]+=3;
	pave->eBins[eIndex+2]+=1;
      }

      break;
    }
  }
  printf("cnt=%d\n",cnt);
}


void loadIdeal(state* sideal,FILE* ifile){
  int N=sideal->N,dummy;
  float x,y,z;
  for(int i=0;i<N;i++){
    dummy=fscanf(ifile,"%f %f %f\n",&x,&y,&z);
    sideal->x[3*i+1]=x;
    sideal->x[3*i+2]=y;
    sideal->x[3*i+3]=z;
  }
}

void resetWindow(std::queue<int>* accepts,float* acceptAvg, int aLen){
  for(int i=0;i<aLen/2;i++){
    accepts->pop();
    accepts->pop();
    accepts->push(1);
    accepts->push(0);
  }
  *acceptAvg=0.5;
}

void printPave(FILE *fp,ARGST* args){
  fprintf(fp,"%f %f %d\n",args->minE,args->delE,args->nBin);
  printf("%d\n",args->nBin);
  for(int i=0;i<args->nBin;i++)
    fprintf(fp,"%d\n",args->eBins[i]);  
}
