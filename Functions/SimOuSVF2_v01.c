/* Change this to reflect the apropriate header file */
#include <math.h>
#include <limits.h>
#include "mex.h"
#include "matrix.h"


double sexp(double x) 
{ 

    return x<log(1.5)?exp(x):1.5*sqrt(1-log(1.5)+pow(x,2)/log(1.5)); 
} 

void SimOuSVF2(double *parm, double *eB1, double *eB2, double *eWp, double *v1, double *v2, double *lpr, double *NSecs, double *Nd)
{
    mwIndex i, nSecs;
	double del, mu, be0, be1, be2, al1, al2, phi, rho1, rho2, scale;
    mu   = parm[0];
	be0  = parm[1];
	be1  = parm[2];
	be2  = parm[3];
    al1  = parm[4];
    al2  = parm[5];
    phi  = parm[6];
    rho1 = parm[7];
    rho2 = parm[8];
    del = 1/Nd[0];
	nSecs = (mwIndex)NSecs[0];
    
	for (i=0; i<(nSecs-1); i++) {

        scale      = del*sexp(be0+be1*v1[i]+be2*v2[i]);	
        
        lpr[i+1]   = lpr[i] + mu*del + sqrt(scale)*eWp[i];
        
        v1[i+1]    = exp(al1*del)*v1[i] + sqrt((1-exp(2*al1*del))/(-2*al1))*eB1[i];

        v2[i+1]    = (1+al2*del)*v2[i] + (1+phi*v2[i])*sqrt(del)*eB2[i];
 	}
}

/* The gateway routine */
void mexFunction( int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[])
{
	/*[lpr v] = SimOuSV_BNS(NSecs, Nd, parameters, v0, lpr0, e)*/

    double delta;
	double *NSecs, *Nd, *parameters, *v0, *lpr0, *e, *lpr, *v1, *v2, *eB1, *eB2, *eWp;

    /*  Check for proper number of arguments. */
	if(nrhs!=8)
		mexErrMsgTxt("Six inputs required.");
	if(nlhs>3)
		mexErrMsgTxt("One or two outputs only.");

	/*  Create a pointer to the input matrices . */
	NSecs      = mxGetPr(prhs[0]);
	Nd         = mxGetPr(prhs[1]);
    parameters = mxGetPr(prhs[2]);
    v0         = mxGetPr(prhs[3]);
    lpr0       = mxGetPr(prhs[4]);
    eB1        = mxGetPr(prhs[5]);
    eB2        = mxGetPr(prhs[6]);
    eWp        = mxGetPr(prhs[7]);
	
	/*  Set the output pointer to the output matrix. */
	plhs[0] = mxCreateDoubleMatrix((int)NSecs[0],1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix((int)NSecs[0],1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix((int)NSecs[0],1, mxREAL);
	

	/*  Create a C pointer to a copy of the output matrix. */
	lpr = mxGetPr(plhs[0]);
    v1  = mxGetPr(plhs[1]);
    v2  = mxGetPr(plhs[2]);
	
    v1[0]  = (mxGetM(prhs[3])!=0?v0[0]:1);
    v2[0]  = (mxGetM(prhs[3])!=0?v0[1]:1);
    lpr[0] = (mxGetM(prhs[4])!=0?lpr0[0]:0.0);

    /*  Call the C subroutine. */
	SimOuSVF2(parameters, eB1, eB2, eWp, v1, v2, lpr, NSecs, Nd);
}

