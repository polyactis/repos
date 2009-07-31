/*=================================================================
 * MexDiagOfTriXTri.c
 *    Implements mex call to DiagOfTriXTri
 *
 *   GADA -- Genome Alteration Detection Algorithm 
 *   Copyright (C) 2008,2009  Childrens Hospital of Los Angeles
 *   Author: Roger Pique-Regi    rpique@gmail.com

  *=================================================================*/

#include "BaseGADA.h"
#include "matlabdefines.h"

void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    
    double *ll; //2D Array containing the [tu';t0';tl';b']
    double *l0;
    double *lu,*rl,*r0,*ru;
    double *d;
    
    int N;     //Number of variables
    
/* Examine input (right-hand-side) arguments. */
    if (nrhs != 6) { 
    	mexErrMsgTxt("Only six input arguments are required.\n"); 
    } 
    if (nlhs != 1) {
    	mexErrMsgTxt("One output argument is required"); 
    } 
    
    N = mxGetN(prhs[1]); 
    ll = mxGetPr(prhs[0]);
    l0 = mxGetPr(prhs[1]);
    lu = mxGetPr(prhs[2]);
    rl = mxGetPr(prhs[3]);
    r0 = mxGetPr(prhs[4]);
    ru = mxGetPr(prhs[5]);
    
    
    /* Create output variables*/
    //mexPrintf("\nThere are %d left-hand-side argument(s).\n", nlhs);
    
  	plhs[0]=mxCreateDoubleMatrix(1,N,mxREAL); //x
    d=mxGetPr(plhs[0]);
    
    DiagOfTridXTrid(ll,l0,lu,rl,r0,ru,d,N);
     // d=DiagOfTriXTri(ll,l0,lu,rl,r0,ru)
 }
