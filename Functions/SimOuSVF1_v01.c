/* Change this to reflect the apropriate header file */
#include <math.h>
#include <limits.h>
#include "mex.h"
#include "matrix.h"

/*
* egarch_core.c -
* This is a helper function and is part of the UCSD_GARCH toolbox
* You can compile it and should work on any platform.
*
* Author: Kevin Sheppard
* kevin.sheppard@economics.ox.ac.uk
* Revision: 3    Date: 3/1/2005
*/
void SimOuSV_BNS_core(double *parameters, double *e, double *v, double *lpr, double *NSecs, double *Nd)
{
    mwIndex i, nSecs;
	double del, mu, be0, be1, alp, rho, rhoConst, scale;
    mu  = parameters[0];
	be0 = parameters[1];
	be1 = parameters[2];
    alp = parameters[3];
    rho = parameters[4];
    del = 1/Nd[0];
    rhoConst = sqrt(1 - rho * rho);
	nSecs = (mwIndex)NSecs[0];
    
	for (i=0; i<(nSecs-1); i++) {

        e[i+nSecs] = rho*e[i] + rhoConst * e[i+nSecs];
        scale      = del*exp(2*(be0+be1*v[i]));	

        lpr[i+1]   = lpr[i] + mu*del + sqrt(scale)*e[i+nSecs];
        v[i+1]     = exp(alp*del)*v[i] + sqrt((1-exp(2*alp*del))/(-2*alp))*e[i];
 	}
}

/* The gateway routine */
void mexFunction( int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[])
{
	/*[lpr v] = SimOuSV_BNS(NSecs, Nd, parameters, v0, lpr0, e)*/

    double delta;
	double *NSecs, *Nd, *parameters, *v0, *lpr0, *e, *lpr, *v;

    /*  Check for proper number of arguments. */
	if(nrhs!=6)
		mexErrMsgTxt("Six inputs required.");
	if(nlhs>2)
		mexErrMsgTxt("One or two outputs only.");

	/*  Create a pointer to the input matrices . */
	NSecs      = mxGetPr(prhs[0]);
	Nd         = mxGetPr(prhs[1]);
    parameters = mxGetPr(prhs[2]);
    v0         = mxGetPr(prhs[3]);
    lpr0       = mxGetPr(prhs[4]);
    e          = mxGetPr(prhs[5]);
	
	/*  Set the output pointer to the output matrix. */
	plhs[0] = mxCreateDoubleMatrix((int)NSecs[0],1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix((int)NSecs[0],1, mxREAL);
	

	/*  Create a C pointer to a copy of the output matrix. */
	lpr = mxGetPr(plhs[0]);
    v   = mxGetPr(plhs[1]);
	
   
    v[0]   = (mxGetM(prhs[3])!=0?v0[0]:1);
    lpr[0] = (mxGetM(prhs[4])!=0?lpr0[0]:0.0);

    /*  Call the C subroutine. */
	SimOuSV_BNS_core(parameters, e, v, lpr, NSecs, Nd);
}

