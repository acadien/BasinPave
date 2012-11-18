#include <math.h>
#define NRANSI
#include "nrutil.h"
#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY2 1.0e-20
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

void mnbrak(float *ax, float *bx, float *cx, float *fa, float *fb, float *fc,
	    float (*func)(float,void*),void* args)
{
	float ulim,u,r,q,fu,dum;

	*fa=(*func)(*ax,args);
	*fb=(*func)(*bx,args);
	if (*fb > *fa) {
		SHFT(dum,*ax,*bx,dum)
		SHFT(dum,*fb,*fa,dum)
	}
	*cx=(*bx)+GOLD*(*bx-*ax);
	*fc=(*func)(*cx,args);
	while (*fb > *fc) {
		r=(*bx-*ax)*(*fb-*fc);
		q=(*bx-*cx)*(*fb-*fa);
		u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
			(2.0*SIGN(FMAX(fabs(q-r),TINY2),q-r));
		ulim=(*bx)+GLIMIT*(*cx-*bx);
		if ((*bx-u)*(u-*cx) > 0.0) {
		  fu=(*func)(u,args);
			if (fu < *fc) {
				*ax=(*bx);
				*bx=u;
				*fa=(*fb);
				*fb=fu;
				return;
			} else if (fu > *fb) {
				*cx=u;
				*fc=fu;
				return;
			}
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func)(u,args);
		} else if ((*cx-u)*(u-ulim) > 0.0) {
		  fu=(*func)(u,args);
			if (fu < *fc) {
				SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
				  SHFT(*fb,*fc,fu,(*func)(u,args))
			}
		} else if ((u-ulim)*(ulim-*cx) >= 0.0) {
			u=ulim;
			fu=(*func)(u,args);
		} else {
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func)(u,args);
		}
		SHFT(*ax,*bx,*cx,u)
		SHFT(*fa,*fb,*fc,fu)
	}
}
#undef GOLD
#undef GLIMIT
#undef TINY2
#undef SHFT
#undef NRANSI
