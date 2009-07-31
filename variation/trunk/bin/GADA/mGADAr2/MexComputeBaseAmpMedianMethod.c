/*=================================================================
 * MexComputeBaseAmpMedianMethod.c 
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
    int *SegLen;
    int K;
    double *SegAmp,*aux;
    double BaseAmp;
    
    // Checking number of input/output parameters
    if (nrhs != 2){ 
    	mexErrMsgTxt("2 input arguments are required"); 
    }
    
    if (nlhs != 1) {
    	mexErrMsgTxt("1 output arguments are required"); 
    }
       
    /* Getting Input parameters */
    SegLen=(int*)mxGetPr(prhs[0]);    
	SegAmp=mxGetPr(prhs[1]);
    K = mxGetNumberOfElements(prhs[0])-1; //Check SegLen and SegAmp 
   
	BaseAmp=CompBaseAmpMedianMethod(SegLen,SegAmp,K);

    // Preparing output results   
    plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL); //
    aux=mxGetPr(plhs[0]);
	aux[0]=BaseAmp;

            
}


