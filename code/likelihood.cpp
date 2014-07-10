//
//  likelihood.cpp
//  phloxsim
//
//  Created by Robin Hopkins on 3/7/13.
//  Copyright (c) 2013 Robin Hopkins. All rights reserved.
//
//#define SIMDAT
//#define FILESDBGR
//#define LIKEDBGR
#define ZeroLike -100000000000

#include "likelihood.h"

#include <iostream>
using std::cout;
#include <string>
using std::string;
#include <fstream>
using std::ifstream;
#include <sstream>
using std::stringstream;
using std::istringstream;


//likelihood function
double MakeLikelihood (vector<double>& predict,vector<int> obs){
    double LH=0;
    double templike=0;
    for (int i=0; i<obs.size();++i){
        if(predict[i]!=0)
            templike= obs[i]*log(predict[i]); //pow(predict[i],obs[i]);
        else if(predict[i]==0 && obs[i]==0)
            templike=0.0;
        else if(predict[i]==0 && obs[i] !=0)
            templike=ZeroLike;
        else cout<< "Likelihood calculation failure";
        LDBG("Pred "<<predict[i]<<" and Obs "<<obs[i]<<" subliklihood: "<< templike<<'\n');
        SIMDATA(predict[i]<<" ");
        LH+=templike;
    };
    SIMDATA('\n');
    return LH;
    
};

// Makes a DemeOut (simulated results) for a Deme that would align with the observed populations - This interpolates given the two neighboring demes that are actually included in the simulation
DemeOut calcBar(vector<DemeOut>& simDemes, double obsPos, Params* par){
    double left;
    
    // The remainder (decimal) after fraction of position is multiplied to number of Demes in the simulation
    double gain=modf((obsPos+par->startPop )/ par->deltaDeme , &left);
    //if(gain < 0) { std::cerr<<"Gain < zero. check startdeme distance. Exiting.\n"; exit(1);}
    LDBG("gain::"<<gain<<" obspos "<<obsPos<<" ");
    
    int leftD=left;                             // Simulated deme left of the observed population
    int rightD=leftD+1;                         // Simulated deme right of the observed populations
    vector<double> phenoBar;
    vector<vector<double> > kBar;
    kBar.resize(simDemes[leftD].ktable.size());
    
    //Interpolate the phenotypes (from the phenoBar) and the values in the K matrix
    for (int i=0; i<simDemes[leftD].simPheno.size(); ++i){
        double pfreq=interpolate(simDemes[leftD].simPheno[i], simDemes[rightD].simPheno[i], gain, 1.0);
        LDBG("left deme "<<leftD<<" freq: "<<simDemes[leftD].simPheno[i]<<" right d "<<rightD<<" freq "<<simDemes[rightD].simPheno[i]<<'\n');
        phenoBar.push_back(pfreq);
        LDBG("interpolated freq: "<<pfreq<<"  ");
        for (int j=0; j<simDemes[leftD].ktable[i].size(); ++j){
            double kfreq=interpolate(simDemes[leftD].ktable[i][j], simDemes[rightD].ktable[i][j], gain, 1.0);
            kBar[i].push_back(kfreq);
        }
    }
    LDBG('\n');
    DemeOut out(phenoBar, kBar);
    return out;                         // This is the interpolated DemeOut that will be used to calc likelihood
    
};

vector<ObsData> readObsFiles(int dataCode){
    
    vector<ObsData> out;
    std::map<int, int> popNames;
    int popNum=0;
    
    stringstream inPops;
    inPops<<"inPhlox_Phen_T"<<dataCode<<".data";
    ifstream popfile(inPops.str().c_str());
    if (!popfile.is_open()) cout<<"Loading "<<inPops.str().c_str()<<" failed\n";
    
    string Loader1;
    getline(popfile,Loader1);
    FILESDBG("\nColumn names in popFile "<<Loader1<<'\n');
    
    while (getline(popfile,Loader1)) {
        ObsData currentPop;
        
        istringstream dataLine(Loader1);
        int pop_Code; dataLine >> pop_Code;
        FILESDBG("Read population "<<pop_Code<<", given PopNum "<<popNum<<": ");
        popNames[pop_Code]=popNum;
        ++popNum;
        
        dataLine >>currentPop.position;
        FILESDBG("Position "<<currentPop.position<<", and PhenoCounts -> ");
        
        int intPlace;
        while(dataLine >> intPlace){
            currentPop.PhenoCounts.push_back(intPlace);
            FILESDBG(intPlace<<" ");
        }
        out.push_back(currentPop);
        FILESDBG('\n');
    }
    
    stringstream inMoms;
    inMoms<<"inPhlox_Fam_T"<<dataCode<<".data";
    ifstream momfile(inMoms.str().c_str());
    if (!momfile.is_open()) cout<<"Loading "<<inMoms.str().c_str()<<" failed\n";
    
    string Loader2;
    getline(momfile,Loader2);
    FILESDBG("\nColumn names in momFile "<<Loader2<<'\n');
    
    while (getline(momfile,Loader2)) {
        FamData myFam;
        istringstream dataLine(Loader2);
        int pop; dataLine >> pop;
        dataLine >> myFam.momPheno;
        
        FILESDBG("Read Fam from Pop "<<pop<<", momPheno "<< myFam.momPheno <<", OffsCounts -> ");
        
        int intPlace;
        while(dataLine >> intPlace){
            myFam.offCount.push_back(intPlace);
            FILESDBG(intPlace<<" ");
        }
        out[popNames[pop]].momOff.push_back(myFam);
        FILESDBG(", went into PopNum "<<popNames[pop]<<'\n');
    }
    momfile.close();
    return out;
};

// Calculates total likelihood
double calcTotLike (vector<DemeOut>& simOut, vector<ObsData>& transObs, Params* par){
    SIMDATA('\n');
    double LH1=0.0;
    double LH2=0.0;
    for (int i=0; i<transObs.size(); ++i){                         // loop through each observed population
        SIMDATA(i<<" ");
        DemeOut myDemeBar=calcBar(simOut,transObs[i].position, par);    // Deside what the "expected data" for sim deme is
        
        //first part of likelihood equation just vectors of expected phenotype frequency and observed phenotype frequency
        LDBG("LH1 for pop "<<i<<"::\n");
        LH1 += MakeLikelihood (myDemeBar.simPheno, transObs[i].PhenoCounts);
        
        for(int j=0; j<transObs[i].momOff.size(); ++j){            //loop through all moms in a population
            int mPheno=transObs[i].momOff[j].momPheno;             // Mom's phenotype for mom j in pop i
            //second part of likelihood equation vector of expected frequency for a mom phenotype and counts from moms
            LDBG("LH2::\n");
           // LH2 += MakeLikelihood(myDemeBar.ktable[mPheno], transObs[i].momOff[j].offCount);
        };
        
    };
    
    SIMDATA('\n');
    for (int i=0; i<transObs.size(); ++i){
        DemeOut myDemeBar=calcBar(simOut,transObs[i].position, par);    // Deside what the "expected data" for sim deme is
        for(int j=0; j<transObs[i].momOff.size(); ++j){            //loop through all moms in a population
            int mPheno=transObs[i].momOff[j].momPheno;             // Mom's phenotype for mom j in pop i
            //second part of likelihood equation vector of expected frequency for a mom phenotype and counts from moms
            LDBG("LH2::\n");
            SIMDATA(i <<" "<<mPheno<<" ");
            LH2 += MakeLikelihood(myDemeBar.ktable[mPheno], transObs[i].momOff[j].offCount);
        };
    };
    //std::cout <<"lh1 "<<LH1<<" lh2 "<<LH2<<'\n';
    if(LH1+LH2 > 0) {
        //std::cerr<<"Positive Log(L) is suspicious. Exiting.\n"; exit(1)
        ;}
    return LH1+LH2;         // total likelihood across all moms and all phenotypes in all populations.
    
       
};
