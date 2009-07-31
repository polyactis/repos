/*=================================================================
 * MexIextWextToSegAmp.c   
 *
 *  SegAmp = MexIextWextToSegAmp( uint32(Iext) , Wext );
 *  SegAmp=IextWextToSegAmp(Iext,Wext);
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
    double *SegAmp;
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
   
    // Preparing output results   
    plhs[0]=mxCreateDoubleMatrix(1,K+1,mxREAL); //
    SegAmp= mxGetPr(plhs[0]);


//    plhs[i]= mxCreateNumericMatrix(1,K+1,mxUINT32_CLASS,mxREAL); //I
//    SegLen= (int*) mxGetPr(plhs[i++]);

	IextWextToSegAmp(Iext,Wext,SegAmp,K);
    


            
}


