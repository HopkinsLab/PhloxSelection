//  This is the header file for the WORLD class which is the main "player" in the simulation (everything is created and exists in a world).

//  World.h
//  phloxsim
//
//  Created by Robin Hopkins on 10/23/12.
//  Copyright (c) 2012 Robin Hopkins. All rights reserved.
//

#ifndef __phloxsim__World__
#define __phloxsim__World__

// These are our general includes from the stl.
#include <iostream>
#include <vector>
using std::vector;

// These are our program specific includes 
#include "Deme.h"
#include "likelihood.h"
#include "Params.h"
#include "utilities.h"
#include <boost/shared_ptr.hpp>
using boost::shared_ptr;


//#define STARTdEME 1.5        //The distance from the start of the simulated transect to the first observed population
//#define DELTAdEME .15        //The distance between each deme in kilometers

// Simulation results for a given deme

//This struct formats the input data. Takes the observed phenotype counts from the pops and turns it into deme position and phenotype frequencies
struct InputPop{
    double position;
    double posDeme;
    vector<int> PhenoCounts;
    vector<double> PhenoFreq;
    
    
    InputPop(double pos, vector<int> phenC, Params* par){
        position=pos;                           //position in KM
        PhenoCounts=phenC;                      //phenotype counts (DB, DR, LB, LR)
        posDeme=(position+ par->startPop )/ par->deltaDeme - 1; //calculates position amoung transect of demes that will be simulated
        vector<double> temp;
        temp.resize(phenC.size());
        double totphen=0;
        //Calculates phenotype frequencies from phenotype counts
        for (int i=0; i<phenC.size(); ++i){
            totphen+=phenC[i];
        };
        
        //std::cout<<"totpheno: "<<totphen<<" ";
        
        for (int j=0; j<phenC.size(); ++j){
            double freq;
            freq=(phenC[j]/totphen);
            temp[j]=freq;
            
        };
        PhenoFreq=temp;
        //std::cout<<"freq: "<<temp[0]<<" "<<temp[1]<<" "<<temp[2]<<" "<<temp[3]<<" "<<'\n';
    };
};

class World{
private:
    vector<Deme> transect;                          // The list of Demes (populations)
    vector<vector<ParentPair> > parents;            // For each offspring genotype there is a list of possible parents
    vector<vector< double> > wX;                    // For each selected region (zone) there is a list of selection coef. (In other parts of the code we assume 2 zones)
    vector<double> pX;                              // Innitial genotype frequencies (from startP out of intitialfreq.cpp)
    vector<vector<double> > startP;                 // Initial genotype freqs for all demes
    gps_type gps;                                   // The Genotype map (see genomap.h)
    double a;                                       // The fraction of assortative mating - same for all demes and based on phenotype.
    unsigned long nGen;
    
public:
    World(Params* par, int key);                         // Creates the world and will exicute the following functions:
    void makeGPS();
    void build_wX(Params* par);
    void buildParents();
    void buildTransect(Params* par, int key);
    void buildNeighborhood(double& m, double&delta, int&migDis);          // For each deme, this identifies neighboring demes that contribute migrants
    
    // These are the functions for setting up the initial genotype freq conditions (initialfreq.cpp)
    vector<double> phenoGenoHW(vector<double> phenoFreq); // Turns phenotype freq into genotype based on HWE
    void initialFreqs(vector<InputPop> obsPop, int totDemes);        // calculates phenotype freqs for all demes
    vector<InputPop> readinput(Params* par, int key);               // Inputs the observed pops for initial conditions
    
    // After building these are the functions executed by simulate across each deme -These functions control the lifecycle and tell each deme what to do (don't actually do the function implied)
    void simulate(Params* paramVals, int key);
    void selectW();
    void migrateMalesW();           
    void matingW();
    double getMaxChg();
    void makeKW();

    // After running simulation these are the output functions (WorldOutput.cpp)
    void outputP();
    void outputPS();
    void outputPM();
    void outputFreq();
    void cout_Gens();
    double output_Gens();
    void outputMaxDelta();
    vector<vector<double> > outputPhenoF();
    vector<vector<vector< double> > >  outputK();
    vector<DemeOut> MakesimDemeOut();
    bool contSim(vector<ObsData>& transObs, Params* par);
};
#ifdef PARAMDBGR
#define PARDBG(x) std::cout << x
#else
#define PARDBG(x)
#endif

#ifdef FILESDBGR
#define FILESDBG(x) std::cout << x
#else
#define FILESDBG(x)
#endif

#ifdef LIKEDBGR
#define LDBG(x) std::cout << x
#else
#define LDBG(x)

#ifdef MIGDBGR
#define MIGDBG(x) std::cout << x
#else
#define MIGDBG(x)
#endif

#endif
#endif /* defined(__phloxsim__World__) */
