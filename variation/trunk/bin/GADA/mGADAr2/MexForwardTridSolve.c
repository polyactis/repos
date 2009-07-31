/*=================================================================
 * MexForwardTridSolve.c 
 *   Implements the mex wrapper for ForwardTridSolve // BackwardTridSolve
 *    [Fx,Fu,Fd,Fl] = ForwardTridSolve(tu,tc,tl,b);
 *    [Bx,Bu,Bd,Bl] = MexBackwardTridSolve(tu,tc,tl,b);  
 *  This is how it will be called from Matlab
 * see also ForwardTridSolve.m
 *
 *
function [x,u,d,l] = ForwardTridSolve(tu,tc,tl,b);
% ForwardTridSolve - Solves Tridiagonal system Tx=b by LDU factorization
% The LDU factorization is obtained by Forward Gaussian elimination
% withouth pivoting
%
% Syntax:  
%   [x,u,d,l] = ForwardTridSolve(tu,tc,tl,b);
%
% Inputs:
%   tu - Upper diagonal of the Tridiagonal matrix T   
%   tc - Main  diagonal of the Tridiagonal matrix T   
%   tl - Lower diagonal of the Tridiagonal matrix T   
%   b  - Right-side of the equation to solve
% 
% Outputs:
%   u  - Upper diagonal of U in T=LDU
%   d  - Main  diagonal of D in T=LDU
%   u  - Lower diagonal of L in T=LDU
%   x  - Solution of the tridiagonal system.
%
% Example:
%   N=100;
%   tc=randn(1,N);
%   tu=randn(1,N-1); 
%   x=randn(1,N);
%   b=TriSymGaxpy(tc,tu,x)
%   [Fx,Fu,Fd,Fl] = ForwardTridSolve(tu,tc,tu,b);
%   [Bx,Bu,Bd,Bl] = MexBackwardTridSolve(tu,tc,tu,b);
%   Fe=sum(abs(Fx-x))
%   Be=sum(abs(Bx-x))
%   FB=sum(abs(Fx-Bx))
%   FBe=sum(abs(0.5*Fx+0.5*Bx-x))
%   
% Implemented in mex built-in function ForwardTridSolve.c
%  
% See also: MexBackwardTridSolve 

% Author: Roger Pique-Regi
% Ming Hsieh Department of Electrical Engineering, University of Southern California 
% Department of Hemathology/Oncology, Childrens Hospital Los Angeles 
% email: rpique@gmail.com
% Website: http://sipi.usc.edu/~piquereg
% March 2008; 

N=length(b);
assert(length(tc)==N);
assert(length(tu)==N-1);
assert(length(tl)==N-1);

[x,u,d,l] = builtin('ForwardTridSolve',tu,tc,tl,b);
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
//    UDLofTrid(tu,tc,tl,N,u,d,l);
//    TridSolveUsingUDL(u,d,l,N,b,x);
//	ForwardTridSolve(tu,tc,tl,N,u,d,l,b,x);
    LDUofTrid(tu,tc,tl,N,u,d,l);
    TridSolveUsingLDU(u,d,l,N,b,x);
    
}
