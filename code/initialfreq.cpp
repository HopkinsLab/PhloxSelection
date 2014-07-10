//
//  initialfreq.cpp
//  phloxsim
//
//  Created by Robin Hopkins on 3/12/13.
//  Copyright (c) 2013 Robin Hopkins. All rights reserved.
//
//#define FILESDBGR

#include <iostream>
using std::cout;
#include <string>
using std::string;
#include <fstream>
using std::ifstream;
#include <sstream>
using std::stringstream;
using std::istringstream;
#include<math.h>

#include "World.h"

// Calculates genotype from phenotype assuming HardyWeinberg
vector<double> World::phenoGenoHW(vector<double> phenoFreq){
    vector<double>pHW;                                  
    double l;
    double r;
    l=sqrt(phenoFreq[2]+phenoFreq[3]);                      // Light allele freq. from pheno
    r=sqrt(phenoFreq[1]+phenoFreq[3]);                      // Red allele freq. from pheno
    FILESDBG("\n pheno0 "<<phenoFreq[0]<<" pheno1 "<<phenoFreq[1]<<" pheno2 "<<phenoFreq[2]<<" pheno3 "<<phenoFreq[3]<<" \nl freq "<<l<<" r freq "<<r<<'\n');
    for (gps_i g=gps.begin(); g!=gps.end(); ++g) {
        Genotype dummyGeno = (*g).first;
        double initialP;
        double a, b, c, d;
        a=b=c=d=0;
        // Place holders for allele frequencies (a is Mi, b is Mh, c is Fi, d is Fh)
        // If allele at locus is recessive the frequencies is l or r and if dominant its 1-l, 1-r.
        if (dummyGeno.M.In ==y) a=l;
        else if (dummyGeno.M.In ==Y) a=(1-l);
        else {std::cout<<" Check startFreqsHW function:: undefined initial genotype\n";}
        
        if (dummyGeno.M.Hu ==h) b=r;
        else if (dummyGeno.M.Hu ==H) b=(1-r);
        else {std::cout<<" Check startFreqsHW function:: undefined initial genotype\n";}
        
        if (dummyGeno.F.In ==y) c=l;
        else if (dummyGeno.F.In ==Y) c=(1-l);
        else {std::cout<<" Check startFreqsHW function:: undefined initial genotype\n";}
        
        if (dummyGeno.F.Hu ==h) d=r;
        else if (dummyGeno.F.Hu ==H) d=(1-r);
        else {std::cout<<" Check startFreqsHW function:: undefined initial genotype\n";}
        
        initialP=a*b*c*d;
        // Genotype frequency is the product of all the allele frequencies making up genotype
        pHW.push_back(initialP);
    };
    return pHW;
    
};


// read in phenotype frequencies from observed data
vector<InputPop> World::readinput(Params* par, int key){
    vector<InputPop> out;
    std::map<int, int> popNames;
    int popNum=0;
    
    stringstream inPops;
    inPops<<"init_Phlox_IF"<<par->wName[key]<<".txt";
    ifstream popfile(inPops.str().c_str());
    if (!popfile.is_open()) cout<<"Loading '"<<inPops.str().c_str()<<"' failed\n";
    
    string Loader1;
    getline(popfile,Loader1);
    FILESDBG("\nFROM World::readinput()\nColumn names in popFile "<<Loader1<<'\n');
    
    while (getline(popfile,Loader1)){
        double tempPos;
        vector<int> tempPhen;
        
        istringstream dataLine(Loader1);
        int pop_Code; dataLine >> pop_Code;
        FILESDBG("Read population "<<pop_Code<<", given PopNum "<<popNum<<": ");
        popNames[pop_Code]=popNum;
        ++popNum;
        
        dataLine >>tempPos;
        FILESDBG("Position "<<tempPos<<", and PhenoCounts -> ");
        
        int intPlace;
        while(dataLine >> intPlace){
            tempPhen.push_back(intPlace);
            FILESDBG(intPlace<<" ");
            
        };
        InputPop currentDeme (tempPos,tempPhen, par);
        out.push_back(currentDeme);
        FILESDBG('\n');
    };
    return out;
};

//calcaulte intitial frequencies for each deme in the world
void World::initialFreqs(vector<InputPop> obsPop, int totDemes){
    FILESDBG("FROM World::initialFreqs()\n");
    int nbPops=obsPop.size();
    int demeCount=0.0;
    //past the end of the Pops for IF, we assume all demes have the same freqs as the last observed Pop
    double FirstPop=obsPop[0].posDeme;
    FILESDBG("FirstPop "<<FirstPop);
    for(int leftEdge= 0; leftEdge < FirstPop; ++leftEdge){
        vector<double>tempP;
        tempP=phenoGenoHW(obsPop[0].PhenoFreq);
        startP.push_back(tempP);
        demeCount+=1;
        FILESDBG("At beginning of IF: "<<obsPop[0].PhenoFreq[0]<<'\n');
    }
    
    for (int pop=0; pop<nbPops-1; ++pop){                 // Loop through each observed population
        double leftPop=obsPop[pop].posDeme;
        double rightPop=obsPop[pop+1].posDeme;
            
        // Loop through each deme between observed pop and next pop
        for (int i=ceil(leftPop); i<rightPop; ++i){
            vector<double> tempFreq;
            vector<double>tempP;
            tempFreq.resize(obsPop[pop].PhenoFreq.size());
            double gain=i-leftPop;
            double deltapop=rightPop-leftPop;
            
            FILESDBG("Deme "<<i<<" left pop "<< leftPop<<" right Pop "<< rightPop<<" has phenotype frequencies: " );
            
            // Interpolate each phenotype frequencies
            for (int j=0; j<tempFreq.size(); ++j){
                double leftfreq=obsPop[pop].PhenoFreq[j];
                double rghtfreq=obsPop[pop+1].PhenoFreq[j];
                double interFreq=interpolate(leftfreq,rghtfreq, gain, deltapop);
                tempFreq[j]=interFreq;
                
                FILESDBG(tempFreq[j]<<" ");
            };
            tempP=phenoGenoHW(tempFreq);
            startP.push_back(tempP);
            demeCount+=1;
            FILESDBG('\n');
        };
    };
    
    //past the end of the Pops for IF, we assume all demes have the same freqs as the last observed Pop
    double lastPop=obsPop[nbPops-1].posDeme;
    FILESDBG("lastPop "<<lastPop);
    for(int rightEdge= ceil(lastPop); rightEdge < totDemes; ++rightEdge){
        vector<double>tempP;
        tempP=phenoGenoHW(obsPop[nbPops-1].PhenoFreq);
        startP.push_back(tempP);
        demeCount+=1;
        FILESDBG("At end of IF: "<<obsPop[nbPops-1].PhenoFreq[0]<<'\n');
        
    };
    if(demeCount !=totDemes){
        std::cerr<<"demeCount is wrong check innitialization "<<demeCount<<" exiting";
        exit(1);}

};


