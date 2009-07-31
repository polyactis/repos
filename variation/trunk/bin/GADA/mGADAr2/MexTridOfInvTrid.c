/*=================================================================
 * MexTridOfInvTrid.c 
 *   Mex call of TridOfInvTrid
 *   Matlab calling [itl,it0,itu] = TridOfInvTrid(Fu,Fd,Bl)
 *
 *
 *   GADA -- Genome Alteration Detection Algorithm 
 *   Copyright (C) 2008,2009  Childrens Hospital of Los Angeles
 *   Author: Roger Pique-Regi    rpique@gmail.com

 *
  *=================================================================*/

#include "BaseGADA.h"
#include "matlabdefines.h"

void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    // input parameters to the main function: tl,t0,d,e
    
    double *Fu,*Fd,*Bl; 
    double *it0;
    double *itu;
    double *itl;
    
    
    int N;     //Number of variables
      
    /* Examine input (right-hand-side) arguments. */
    if (nrhs != 3)  
    	mexErrMsgTxt("Only three input arguments are required.\n"); 
    if (nlhs != 3) 
    	mexErrMsgTxt("Three output arguments are required"); 
    
    N = mxGetNumberOfElements(prhs[1]); 
    
    if (mxGetNumberOfElements(prhs[0]) != N-1 || mxGetNumberOfElements(prhs[2]) != N-1){
        mexErrMsgTxt("Array dimensions not correct.\n");
    }
    
    Fu = mxGetPr(prhs[0]);
    Fd = mxGetPr(prhs[1]);
    Bl = mxGetPr(prhs[2]);
    
    /* Create output variables*/
    //mexPrintf("\nThere are %d left-hand-side argument(s).\n", nlhs);
    
  	plhs[0]=mxCreateDoubleMatrix(1,N-1,mxREAL); 
    itl=mxGetPr(plhs[0]);
    
  	plhs[1]=mxCreateDoubleMatrix(1,N,mxREAL); 
    it0=mxGetPr(plhs[1]);
    
  	plhs[2]=mxCreateDoubleMatrix(1,N-1,mxREAL); 
    itu=mxGetPr(plhs[2]);

    TridOfInvTrid(Fu,Fd,Bl,itl,it0,itu,N);
    
}
