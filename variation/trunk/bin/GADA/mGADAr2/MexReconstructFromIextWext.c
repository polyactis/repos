/*=================================================================
 * MexReconstructFromIextWext.c   
 * 
 * xRec=ReconstructFromIextWext(Iext,Wext);
 * void ReconstructFromIextWext(
		int *Iext ///(in) Breakpoint set defining the segmentation.
		,double *Wext ///(in) Breakpoint coefficients defining the segmentation.
		,double *xRec ///(out) Full reconstruction 
		,int K ///(in) Number of breakpoints
		)
 *
 *   GADA -- Genome Alteration Detection Algorithm 
 *   Copyright (C) 2008,2009  Childrens Hospital of Los Angeles
 *   Author: Roger Pique-Regi    rpique@gmail.com

  *=================================================================*/

#include "BaseGADA.h"
#include <math.h>
#include "matlabdefines.h"
//#include "mex.h"
//#iclude "matrix.h"
             
void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{  
    /* VARIABLES DECLARATION */
    int *Iext;
    int K,k;
    double *xRec;
    double *Wext;
   
    // Checking number of input/output parameters
    if (nrhs != 2){ 
    	mexErrMsgTxt("2 input arguments are required"); 
    }
    
    if (nlhs != 1) {
    	mexErrMsgTxt("1 output arguments are required"); 
    }

       
    /* Getting Input parameters */
    Iext=(int*)mxGetPr(prhs[0]);    
	Wext=mxGetPr(prhs[1]);
    K = mxGetNumberOfElements(prhs[0])-2; //Check Iext and Wext   
   
    //mexPrintf("K%d,M%d\n",K,Iext[K+1]);
    
    // Preparing output results   
    plhs[0]=mxCreateDoubleMatrix(1,Iext[K+1],mxREAL); //
    xRec= mxGetPr(plhs[0]);


	ReconstructFromIextWext(Iext,Wext,xRec,K);
            
}


