//
//  Params.h
//  phloxsim


// Struct of parameter values
// Robin & Mark 19-X-12

#ifndef __phloxsim__paramS__
#define __phloxsim__paramS__

#include <iostream>
#include <vector>
using std::vector;



struct Params{
    //Params for worldS    
    int nWorlds;    // Number of worlds
    vector<int> wName; //world name (number corresponding to transect in world)
    
    // Params for the life cycle:
    vector< vector <double> > w;
    // Selection coefficients, each zone has a vector of coefs the order is DB, DR, LR (alphabetical) within a zone
    double m;           // Migration rate
    int migDis;        //  code corresponding to likelihood distribution for migration (1=normal, 2=geometric)
    double a;           // Assortative mating rate
    
    
    // Params describing the demes:
    int totDemes;       // Number of demes
    int nZones;         // Number of zones
    vector <double> border;      //Input the border as a fraction of the total transect
    int dataInt;
    double startPop;   //distance from the start of transect to start of first observed population
    double deltaDeme; //distance between demes in simulation in km
    
    // Simulation parameters:
    int maxGens;        // Max number of generations
    double maxChange;   // Max change for equalibrium
    
    vector<double> borderRange;
    vector<double> migRange;
    //vector<double> aRange;
    vector<double> wRange;
    
    //simulated annealing parameters
    double tolerance;
    double temperature;
    double temp_step;
    double delta;
    
    
    //MCMC parameter
    vector<double> sdparams;
    
    
    Params();
    Params(const char* inputfile);
    void addUntransformedPars(vector<double> parSet);
    vector<double> VecInverseLogit(vector<double> p);
};

#ifdef PDBGR
#define PDBG(x) std::cout << x << std::endl;
#else
#define PDBG(x)
#endif

#endif /* defined(__phloxsim__paramS__) */
