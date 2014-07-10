//
//  likelihood.h
//  phloxsim
//
//  Created by Robin Hopkins on 3/7/13.
//  Copyright (c) 2013 Robin Hopkins. All rights reserved.
//

#ifndef __phloxsim__likelihood__
#define __phloxsim__likelihood__
//These are our general includes form the stl.
#include <iostream>
#include<vector>
#include<math.h>
using std::vector;

//These are our program specific includes
#include "World.h"                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
#include "utilities.h"
#include "Params.h"

//likelihood function
double MakeLikelihood (vector<double>& predict,vector<int> obs);

DemeOut calcBar(vector<DemeOut>& simDemes, double obsPos, Params* par);

double calcTotLike (vector<DemeOut>& simOut,vector<ObsData>& transObs, Params* par);
vector<ObsData> readObsFiles(int dataCode);

#endif /* defined(__phloxsim__likelihood__) */

#ifdef SIMDBGR
#define SIMDBG(x) std::cout << x
#else
#define SIMDBG(x)
#endif

#ifdef SIMDAT
#define SIMDATA(x) std::cout<< x
#else
#define SIMDATA(x)
#endif