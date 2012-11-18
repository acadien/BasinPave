/* note #undef's at end of file */
#define NRANSI
#include "nrutil.h"
#define TOL 2.0e-4

int ncom;
float *pcom,*xicom,(*nrfunc)(float [],void*);

void linmin(float p[], float xi[], int n, float *fret, float (*func)(float [], void*), void* args)
{
	float brent(float ax, float bx, float cx,
		    float (*f)(float,void*), float tol, float *xmin,void* args);
	float f1dim(float x, void* args);
	void mnbrak(float *ax, float *bx, float *cx, float *fa, float *fb,
		    float *fc, float (*func)(float,void*),void* args);
	int j;
	float xx,xmin,fx,fb,fa,bx,ax;

	ncom=n;
	pcom=vector(1,n);
	xicom=vector(1,n);
	nrfunc=func;
	for (j=1;j<=n;j++) {
		pcom[j]=p[j];
		xicom[j]=xi[j];
	}
	ax=0.0;
	xx=1.0;
	mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim,args);
	*fret=brent(ax,xx,bx,f1dim,TOL,&xmin,args);
	for (j=1;j<=n;j++) {
		xi[j] *= xmin;
		p[j] += xi[j];
	}
	free_vector(xicom,1,n);
	free_vector(pcom,1,n);
}
#undef TOL
#undef NRANSI
