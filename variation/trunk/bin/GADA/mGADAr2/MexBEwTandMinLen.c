/*=================================================================
 * MexSBL_BE.c   
 *   Generates the Mex for the SBLandBE function in BaseGADA.c
 * Debugging defines:
 *   _DEBUG_SBLBE_
 *
 *   GADA -- Genome Alteration Detection Algorithm 
 *   Copyright (C) 2008,2009  Childrens Hospital of Los Angeles
 *   Author: Roger Pique-Regi    rpique@gmail.com
  *=================================================================*/

/*
  This File is part of GADA

  GADA v2.0 Genome Alteration Detection Algorithm 
  Copyright (C) 2008,2009  Childrens Hospital of Los Angeles

  GADA is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  any later version.

  GADA is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with GADA.  If not, see <http://www.gnu.org/licenses/>.

  Authors: 
  Roger Pique-Regi    piquereg@usc.edu
*/


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
    int *Iext;
    int K;
    int k;
    int MinSegLen;
    double *Wext;
    int *auxIext;
    double *auxWext;
    
    // Checking number of input/output parameters
    if (nrhs != 5){ 
    	mexErrMsgTxt("5 input arguments are required"); 
    }
    
    if (nlhs != 2) {
    	mexErrMsgTxt("2 output arguments are required"); 
    }

       
    /* Getting Input parameters */
    auxIext=(int*)mxGetPr(prhs[0]);    
    auxWext=mxGetPr(prhs[1]);  
    sigma2=mxGetScalar(prhs[2]);    
    T=mxGetScalar(prhs[3]);
    MinSegLen=(int)mxGetScalar(prhs[4]);  

    K = mxGetNumberOfElements(prhs[0])-2;    
	
	Iext=myCalloc(K+2,sizeof(int));
	Wext=myCalloc(K+1,sizeof(double));

	for(k=0;k<=K+1;++k)
		Iext[k]=auxIext[k];
	for(k=0;k<=K;++k)
		Wext[k]=auxWext[k];

	BEwTandMinLen(Wext,Iext,&K,sigma2,T,MinSegLen);

    // Preparing output results
    plhs[0]= mxCreateNumericMatrix(1,K+2,mxUINT32_CLASS,mxREAL); //I
  	auxIext= (int*) mxGetPr(plhs[0]);
    
    plhs[1]=mxCreateDoubleMatrix(1,K+1,mxREAL); //w
    auxWext= mxGetPr(plhs[1]);
        
	for(k=0;k<=K+1;++k)
          auxIext[k]=Iext[k];
	for(k=0;k<=K;++k)
          auxWext[k]=Wext[k];    
}


