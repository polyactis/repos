/*=================================================================
 * MexCollapseAmpTtest.c  
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
    double sigma2;
    double T;
	double BaseAmp;
    int *SegLen;
    double *SegAmp;
    int L;
    double *OutSegAmp;
    
    // Checking number of input/output parameters
    if (nrhs != 5){ 
    	mexErrMsgTxt("5 input arguments are required"); 
    }
    
    if (nlhs != 1) {
    	mexErrMsgTxt("1 output arguments are required"); 
    }

       
    /* Getting Input parameters */
    SegLen=(int*)mxGetPr(prhs[0]);    
    SegAmp=mxGetPr(prhs[1]); //NotNecessary 
	BaseAmp=mxGetScalar(prhs[2]);  
    sigma2=mxGetScalar(prhs[3]);    
    T=mxGetScalar(prhs[4]);

    L = mxGetNumberOfElements(prhs[0]);    

    // Preparing output results
	plhs[0]=mxDuplicateArray(prhs[1]);
    OutSegAmp= mxGetPr(plhs[0]);

	CollapseAmpTtest(OutSegAmp,SegLen,L,BaseAmp,sigma2,T);
			

}


