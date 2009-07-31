/*=================================================================
 * matlabdefines.h 
 *
 *   GADA -- Genome Alteration Detection Algorithm 
 *   Copyright (C) 2008,2009  Childrens Hospital of Los Angeles
 *   Author: Roger Pique-Regi    rpique@gmail.com
 *=================================================================*/
#ifndef _MATLABDEFINES_H_
#define _MATLABDEFINES_H_

//Matlab printing and memory allocation 
#include "mex.h"
#define myPrintf mexPrintf //Syntax of myPrintf similar to that of fprintf('asfd',var,var)
//memory allocation in Matlab
//#define myDoubleMAlloc(mysize) mxGetPr(mxCreateDoubleMatrix(1,mysize,mxREAL))
//#define myIntMalloc(mysize)    mxGetPr(mxCreateNumericMatrix(1,mysize,mxUINT32_CLASS,mxREAL))
#define myCalloc(n,size) mxCalloc(n, size)
#define myRealloc(ptr,size) mxRealloc(ptr,size) // Checking if this fails or causes memory leak...// checked and it does not seem it is this....
#define myFree(ptr) mxFree(ptr)

#endif //_MATLABDEFINES_H_
