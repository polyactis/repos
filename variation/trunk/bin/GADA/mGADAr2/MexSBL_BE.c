/*=================================================================
 * MexSBL_BE.c   
 *   Generates the Mex for the SBLandBE function in BaseGADA.c
 * [Iext,Wext]=MexSBL_BE(y,sigma2,a,T,MinSegLen);
 * Debugging defines:
 *   _DEBUG_SBLBE_
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
    double a;
    double T;
    int *Iext;
    int K;
    int i=0;
    int M ;
//    int NumEMsteps;    
    int MinSegLen;
//     double *alpha;
    double *Wext;
//     double *sigw;
//     double tol;
//     int maxit;
//     double maxalpha;
    double *y; //array containing noisy signal to find the coefficients
    int *outIext;
    double *outWext;
    
    // Checking number of input/output parameters
    if (nrhs != 5){ 
    	mexErrMsgTxt("5 input arguments are required"); 
    }
    
    if (nlhs != 2) {
    	mexErrMsgTxt("2 output arguments are required"); 
    }

       
    /* Getting Input parameters */
    y=mxGetPr(prhs[0]);    
    sigma2=mxGetScalar(prhs[1]);    
    a=mxGetScalar(prhs[2]);  
    T=mxGetScalar(prhs[3]);
    MinSegLen=(int)mxGetScalar(prhs[4]);
    ////I=(int*)mxGetPr(prhs[4]);
    //tol=mxGetScalar(prhs[5]);
    //maxalpha=mxGetScalar(prhs[6]);
    //maxit=mxGetScalar(prhs[7]);    

    M = mxGetNumberOfElements(prhs[0]);    

    //// K = mxGetNumberOfElements(prhs[4]);    
    //// alpha=mxGetPr(mxCreateDoubleMatrix(1,K,mxREAL));  //
    //// w=mxGetPr(mxCreateDoubleMatrix(1,M,mxREAL));    
    //// sigw=mxGetPr(mxCreateDoubleMatrix(1,M,mxREAL));
    
    #ifdef _DEBUG_SBLBE_
        myPrintf("_SBLBE_ Precal\n");
    #endif
    K=SBLandBE(y,M,&sigma2,a,T,MinSegLen,&Iext,&Wext);
    
    #ifdef _DEBUG_SBLBE_
        myPrintf("_SBLBE_ K = %d \n",K);
    #endif

    // Preparing output results
    i=0;
    plhs[i]= mxCreateNumericMatrix(1,K+2,mxUINT32_CLASS,mxREAL); //I
  	outIext= (int*) mxGetPr(plhs[i++]);
    
    plhs[i]=mxCreateDoubleMatrix(1,K+1,mxREAL); //w
    outWext= mxGetPr(plhs[i++]);
        
    for (i=0;i<K+1;i++){        
          outIext[i]=Iext[i]; 
          outWext[i]=Wext[i];
    }
    outIext[K+1]=Iext[K+1];
    mxAssert(Iext[K+1]==M,"Problem in the Mexfile");

    #ifdef _DEBUG_SBLBE_
       // myPrintf("_SBLBE_ I[0,1,2,K-1,K] %d %d %d %d %d \n",I[0],I[1],I[2],I[K-1],I[K]);
    #endif
    
}


