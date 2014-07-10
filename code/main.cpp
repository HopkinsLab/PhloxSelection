// This file is part of the PHLOX simulation project created by Raf and Ra. This simulation will eventually calculate the likelihood of the observed flower color allele frequencies and ?neutral genetic variation? given parameters for selection, migration, and assortative mating.

//  main.cpp

//  Created by Rafael Guerrero on 10/24/12.
//  Copyright (c) 2012 University of Texas. All rights reserved.
//
// These are includes from standard library that will handle input and output and streams
#include <iostream>
using std::cout;
#include <sstream>
using std::stringstream;
using std::string;

#include "ran.h"
#include <boost/random.hpp>
// Global declarations:

// Random generator recommended by Boost; declared globally so that it can be
//	 seeded in main.cpp and then used elsewhere

using boost::mt19937;
boost::mt19937 gen;
// In real use, seed the random number generator from the system clock -- don't forget to do this!:
//




// These are our project-specific includes
#include "World.h"
#include "likelihood.h"
#include "utilities.h"
#include "mcmc.h"
#define mcmcBurnIn 0
#define mcmcRealz 75000

// This is the main function that inputs parameters, runs simulation and outputs likelihoods or frequencies

int main(int argc, const char * argv[])

{   //gen.seed(2960);
    gen.seed(static_cast<unsigned int>(std::time(0)));
    clock_t init, end;
    init=clock();
     
    Params* GlobalPars=new Params();
    double burnAcceTot = 0;
    double realzAcceTot = 0;
     
    Mchain chain(GlobalPars);
    //chain.mcmcbord(chain.totRuns);
   // chain.mcmcStep();
   // chain.mcmcMultiStep(chain.totRuns);
//    for (int i=0; i<mcmcBurnIn; ++i){
//        burnAcceTot=chain.mcmcStep();
//        chain.dump<<burnAcceTot<<" ";
//        chain.printset();
//        
//    }
//    
    for (int i=0; i<mcmcRealz; ++i){
        realzAcceTot=chain.mcmcStep();
        chain.dump<<realzAcceTot<<" ";
        chain.printset();
    }
    
   return 0;
}

