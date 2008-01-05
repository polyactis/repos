////////////////////////////////////////
///////////MCMCFunction.h///////////////
////////////////////////////////////////

//This function provides a skeleton for MCMC analysis

#ifndef _MCMCFUNCTION_H_
#define _MCMCFUNCTION_H_
/**/
#include <iostream>
#include <math.h>
#include <string>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <istream>
#include <iostream.h>
#include <iomanip>
#include <ctype.h>
#include <time.h>
#include <malloc.h>
#include <memory>
#include <cassert>
#include <stdio.h>

using namespace std;

//extern void ProposeNewState();

//extern double CalculateLnHastingsRatio();

//extern void AcceptNewState();

//extern void RejectNewState();

extern void DoMCMC(int,int);

#endif  //MCMCFUNCTION
