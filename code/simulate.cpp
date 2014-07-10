//
//  simulate.cpp
//  phloxsim
//
//  Created by Robin Hopkins on 10/19/12.
//  Copyright (c) 2012 Robin Hopkins. All rights reserved.
//
//#define SIMDBGR

#define TestGens 5000
#define badLogL -80000
#include "World.h"


//#include </usr/llvm-gcc-4.2/lib/gcc/i686-apple-darwin11/4.2.1/include/omp.h>
//int omp_thread_count() {
//    int n = 0;
//#pragma omp parallel reduction(+:n)
//    n += 1;
//    return n;
//}
//
//void World::simulate(Params* paramVals){
//    nGen = 0;
//    double change = 0.0;
//    do{selectW();                // Do selection in world
//        migrateMalesW();           // Do migration between demes in world
//        matingW();                 // Do mating, meiosis, & fertilization
//        change = getMaxChg();   // Find max genotype freq change across gen
//        nGen++;
//    } while(nGen < paramVals->maxGens && change > paramVals->maxChange);
//        
//        makeKW();
//  
//};

void World::simulate(Params* par, int key){
    nGen = 0;
    
    vector<ObsData> transObs = readObsFiles(par->wName[key] );     // vector of observered data structs from readin file
    
    do{
        selectW();                // Do selection in world
        migrateMalesW();           // Do migration between demes in world
        matingW();                 // Do mating, meiosis, & fertilization
        nGen++;
    } while(contSim(transObs, par));
    //std::cout<<nGen<<'\n';
    makeKW();
    
};

bool World::contSim( vector<ObsData>& transObs, Params* par){
    double change = getMaxChg();
    bool equilibrium = change> par->maxChange;
    bool maxlength = nGen < par->maxGens;
    bool incrLikelihood = true;
    if (nGen % TestGens == 0){
        makeKW();
        vector<DemeOut>testsimOut=MakesimDemeOut();
         double testLike =calcTotLike (testsimOut, transObs, par);
        if (testLike < badLogL)  incrLikelihood=false;
        //std::cout<<"Likelihood after "<<nGen<<" generations: "<<testLike<<'\n';
        };
    return equilibrium && maxlength && incrLikelihood;
};

// Functions to tell each deme to do selection, migration and mating
void World::makeKW(){
    for (int i=0; i!=transect.size(); ++i){
        transect.at(i).makeKD(parents, gps, a);
    }
};
void World::selectW(){
    
//#pragma omp parallel for
    for (int i=0; i<transect.size(); ++i){
       transect[i].selectD();
    }
};

void World::migrateMalesW(){
    for (int i=0; i!=transect.size(); ++i){
        transect.at(i).migrateD();
    }
};

void World::matingW(){
//#pragma omp parallel for
    for (int i=0; i<transect.size(); ++i){
        transect[i].matingD(parents, gps, a);
    }
};

double 	World::getMaxChg(){

double delta_max=0.0;
for (int j =0; j< transect.size(); ++j) {
    delta_max= std::max( transect.at(j).getPDel(), delta_max);
}
return delta_max;
}