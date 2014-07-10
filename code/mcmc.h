//
//  mcmc.h
//  phloxMCMC
//
//  Created by Robin Hopkins on 8/7/13.
//  Copyright (c) 2013 Robin Hopkins. All rights reserved.
//

#ifndef phloxMCMC_mcmc_h
#define phloxMCMC_mcmc_h

#define NPARs 12
#define Mlink 1
#define ParamStart 5


#include <iostream>
using std::cout;
#include <sstream>
using std::stringstream;
using std::string;

#include <boost/random.hpp>

#include "nr3.h"
#include "ran.h"
#include <math.h>
#include <time.h>



#include "likelihood.h"

// struct of parameter values and corresponding likelihood

struct State{
    VecDoub set;
    double likeli;
};

class Mchain{
public:
    int callcounter;
    Params* global;
    vector<vector<ObsData> > dataFromFile;
    std::ofstream dump;
    State current;
    vector<vector<double> > seedvec;        //this will be the list of parameters to run for the MCMC chain
    int totRuns;                            // This is the total number of steps the MCMC chain will take (based on the seeds created by R

    Mchain(Params* GlobalPars);
    double oneRun(VecDoub & pars_of_interest);
    double oneRunPrior(VecDoub & pars_of_interest);
    VecDoub startingSet();
    VecDoub propose(vector<double>& sdparams);
    VecDoub proposeLN(vector<double>& sdpar);
    VecDoub multiPropose(int stepcnt);
    double mcmcStep();
    
    double mcmcMultiStep(int nsteps);   //Call this to start run from already established seedset
    double mcmcbord(int nsteps);
    VecDoub multiBord(int stepcnt);
    void printset();
    vector<vector<double> > inseeds();      // this reads in seeds from R's output and makes the seedvec.

};


#ifdef MCMCDBGR
#define MCMCDBG(x) std::cout << x
#else
#define MCMCDBG(x)
#endif

#ifdef SEEDDBGR
#define SEEDDBG(x) std::cout << x
#else
#define SEEDDBG(x)
#endif



#endif
