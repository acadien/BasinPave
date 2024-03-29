/* note #undef's at end of file */
#define NRANSI
#include "nrutil.h"

extern int ncom;
extern float *pcom,*xicom,(*nrfunc)(float [], void*);

float f1dim(float x, void* args)
{
	int j;
	float f,*xt;

	xt=vector(1,ncom);
	for (j=1;j<=ncom;j++) xt[j]=pcom[j]+x*xicom[j];
	f=(*nrfunc)(xt,args);
	free_vector(xt,1,ncom);
	return f;
}
#undef NRANSI
