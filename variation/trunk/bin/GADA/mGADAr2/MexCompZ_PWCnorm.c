/*=================================================================
 * MexCompZ_PWCnorm.c   
 *
 *  z = MexCompZ_PWCnorm( double(y));
 *    void CompZ(double *y,double *z,int M)
 *
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
    int M;
    double *z;
    double *y;
    
    // Checking number of input/output parameters
    if (nrhs != 1){ 
    	mexErrMsgTxt("1 input arguments are required"); 
    }
    
    if (nlhs != 1) {
    	mexErrMsgTxt("1 output arguments are required"); 
    }
       
    /* Getting Input parameters */
    y = mxGetPr(prhs[0]);    
    M = mxGetNumberOfElements(prhs[0]); 
   
    // Preparing output results   
    plhs[0]=mxCreateDoubleMatrix(1,M-1,mxREAL); //
    z = mxGetPr(plhs[0]);

	CompZ(y,z,M);
            
}


