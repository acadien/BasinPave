#include <math.h>
#define NRANSI
#include "nrutil.h"
#define TINY 1.0e-25
#define ITMAX 1000

void powell(float* p, float **xi, int n, float ftol, int *iter, float *fret,
	    float (*func)(float [], void*),void* args)
{
	void linmin(float p[], float xi[], int n, float *fret,
		    float (*func)(float [], void*),void* args);
	int i,ibig,j;
	float del,fp,fptt,t,*pt,*ptt,*xit;
	//	return;
	pt=vector(1,n);
	ptt=vector(1,n);
	xit=vector(1,n);
	*fret=(*func)(p,args);
	for (j=1;j<=n;j++) pt[j]=p[j];
	for (*iter=1;;++(*iter)) {
		fp=(*fret);
		ibig=0;
		del=0.0;
		for (i=1;i<=n;i++) {
			for (j=1;j<=n;j++) xit[j]=xi[j][i];
			fptt=(*fret);
			linmin(p,xit,n,fret,func,args);
			if (fptt-(*fret) > del) {
				del=fptt-(*fret);
				ibig=i;
			}
		}
		if (2.0*(fp-(*fret)) <= ftol*(fabs(fp)+fabs(*fret))+TINY) {
			free_vector(xit,1,n);
			free_vector(ptt,1,n);
			free_vector(pt,1,n);
			return;
		}
		if (*iter == ITMAX) nrerror("powell exceeding maximum iterations.");
		for (j=1;j<=n;j++) {
			ptt[j]=2.0*p[j]-pt[j];
			xit[j]=p[j]-pt[j];
			pt[j]=p[j];
		}
		fptt=(*func)(ptt,args);
		if (fptt < fp) {
			t=2.0*(fp-2.0*(*fret)+fptt)*SQR(fp-(*fret)-del)-del*SQR(fp-fptt);
			if (t < 0.0) {
			  linmin(p,xit,n,fret,func,args);
				for (j=1;j<=n;j++) {
					xi[j][ibig]=xi[j][n];
					xi[j][n]=xit[j];
				}
			}
		}
	}
}
#undef ITMAX
#undef NRANSI
