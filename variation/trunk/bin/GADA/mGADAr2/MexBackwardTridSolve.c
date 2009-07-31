/*=================================================================
 * MexForwardTridSolve.c 
 *   Implements the mex wrapper for ForwardTridSolve // BackwardTridSolve
 *    [Fx,Fu,Fd,Fl] = MexForwardTridSolve(tu,tc,tl,b);
 *    [Bx,Bu,Bd,Bl] = MexBackwardTridSolve(tu,tc,tl,b);  
 *  This is how it will be called from Matlab
 * see also MexForwardTridSolve.m
 *
 %   GADA -- Genome Alteration Detection Algorithm 
%   Copyright (C) 2008,2009  Childrens Hospital of Los Angeles
%   Author: Roger Pique-Regi    rpique@gmail.com

  *=================================================================*/

#include "BaseGADA.h"
#include "matlabdefines.h"

void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    
    double *tu,*tc,*tl,*b; //2D Array containing the [tu';tc';tl';b']
    double *u,*d,*l,*x;    //Computes the decomposition
    
    int N;     //Number of variables
          
    if (nrhs != 4) { 
    	mexErrMsgTxt("Four input arguments are required.\n"); 
    } 
    if (nlhs != 4) {
    	mexErrMsgTxt("Four output arguments are required"); 
    } 
    
    N = mxGetNumberOfElements(prhs[1]); 
    tu = mxGetPr(prhs[0]); 
    tc = mxGetPr(prhs[1]); 
    tl = mxGetPr(prhs[2]);
    b  = mxGetPr(prhs[3]);
    
  	plhs[0]=mxCreateDoubleMatrix(1,N,mxREAL); //x
    x=mxGetPr(plhs[0]);    
  	plhs[1]=mxCreateDoubleMatrix(1,N-1,mxREAL); //
    u=mxGetPr(plhs[1]);
  	plhs[2]=mxCreateDoubleMatrix(1,N,mxREAL); //
    d=mxGetPr(plhs[2]);
  	plhs[3]=mxCreateDoubleMatrix(1,N-1,mxREAL); //
    l=mxGetPr(plhs[3]);

//    BackwardTridSolve(tu,tc,tl,N,u,d,l,b,x);
    UDLofTrid(tu,tc,tl,N,u,d,l);
    TridSolveUsingUDL(u,d,l,N,b,x);
//	ForwardTridSolve(tu,tc,tl,N,u,d,l,b,x);
//    LDUofTrid(tu,tc,tl,N,u,d,l);
//    TridSolveUsingLDU(u,d,l,N,b,x);
    
}
