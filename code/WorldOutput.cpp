//
//  WorldOutput.cpp
//  phloxsim
//
//  Created by Robin Hopkins on 3/18/13.
//  Copyright (c) 2013 Robin Hopkins. All rights reserved.
//

#include <iomanip>
#include <iostream>
#include <vector>
using std::vector;

#include "World.h"


void World::outputP(){
    
    for(int i=0; i<transect.size();++i){
        vector<shared_ptr<double> > d =transect[i].getP();
        for(int j=0; j<d.size();++j){
            std::cout<<"P: "<<std::setprecision(3)<<*d[j]<<" ";
        }
        std::cout<<'\n';
    }
    
};

// Outputs allelic frequencies
void World::outputFreq(){
    
    for (int i=0; i<transect.size();++i){
        double lfreq=0.0;
        double dfreq=0.0;
        double rfreq=0.0;
        double bfreq=0.0;
        vector<shared_ptr<double> > d =transect[i].getP();
        for(gps_i g=gps.begin(); g!=gps.end(); ++g) {
            Genotype dummyGeno = (*g).first;
            int index=(*g).second;
            
            if (dummyGeno.M.In ==y)lfreq+=((*d[index])*.5);
            else if (dummyGeno.M.In==Y) dfreq+=((*d[index])*.5);
            else {std::cout<<" Check outputFreq for misdefined genotype\n";}
            
            if (dummyGeno.F.In ==y)lfreq+=((*d[index])*.5);
            else if (dummyGeno.F.In==Y) dfreq+=((*d[index])*.5);
            else {std::cout<<" Check outputFreq for misdefined genotype\n";}
            
            if (dummyGeno.M.Hu ==h)rfreq+=((*d[index])*.5);
            else if (dummyGeno.M.Hu==H) bfreq+=((*d[index])*.5);
            else {std::cout<<" Check outputFreq for misdefined genotype\n";}
            
            if (dummyGeno.F.Hu ==h)rfreq+=((*d[index])*.5);
            else if (dummyGeno.F.Hu==H) bfreq+=((*d[index])*.5);
            else {std::cout<<" Check outputFreq for misdefined genotype\n";}
        }
        
        std::cout<<"light frequency "<<lfreq<<" Dark frequency "<<dfreq<<" red frequency "<<rfreq<<" blue frequency "<<bfreq<<'\n';
    }
    
}
// Output genotype freq after selection
void World::outputPS(){
    
    for(int i=0; i<transect.size();++i){
       vector<shared_ptr<double> > d =transect[i].getPS();
        for(int j=0; j<d.size();++j){
            std::cout<<"PS: "<<std::setprecision(3)<<d[j]<<" ";
        }
        std::cout<<'\n';
    }};
// Output genotype freq after migration
void World::outputPM(){
    
    for(int i=0; i<transect.size();++i){
        vector<shared_ptr<double> > d =transect[i].getPM();
        for(int j=0; j<d.size();++j){
            std::cout<<"PM: "<<std::setprecision(3)<<*d[j]<<" ";
        }
        std::cout<<'\n';
    }};

void World::cout_Gens(){
    std::cout<<"number of generations "<<nGen<<'\n';
};

double World::output_Gens(){
    return nGen;
};

void World::outputMaxDelta(){
    std::cout<<"max change in genotype frequency "<<getMaxChg()<<'\n';
};

// Output phenotype frequencies
vector<vector<double> > World::outputPhenoF(){
    vector<vector<double> > simPhenoVec;
    simPhenoVec.resize(transect.size());
    vector<double> pheno;
    pheno.resize(NPHENOS);
    for(int i=0; i<transect.size();++i){
        vector<double> pheno=transect[i].getPheno(gps);
        std::cout<<"Phenotype freq: DB "<< pheno[0]<<" DR "<< pheno[1]<<" LB "<< pheno[2]<<" LR "<< pheno[3]<<'\n';
        simPhenoVec[i]=pheno;
    }
    return simPhenoVec;                             // Vector of expected phenotypes for each deme insimulation
    
};


vector<vector<vector< double> > > World::outputK(){
    vector<vector<vector< double> > >  kTableVec;
    kTableVec.resize(transect.size());
    
    for (int i=0; i<transect.size(); ++i){                      //Loop through all demes in world
        vector<vector<double> > k = transect[i].getK();
        kTableVec[i]=k;
        std::cout<<"population "<< i <<'\n';
        for (int m=0; m<k.size();++m){
            std::cout<<"female with phenotype "<< m <<'\n';
            for (int of=0; of<k[m].size();++of){
                std::cout<<"offspring "<<of<<" probability "<<k[m][of]<<'\n';
                
            };
        };
    };
    return kTableVec;                                // Vector of k tables for each deme in simulation
};

vector<DemeOut> World::MakesimDemeOut(){
    
    DemeOut tempDeme;
    vector<DemeOut> simResults;
    simResults.resize(transect.size());
    
    for(int i=0; i<transect.size(); ++i){
        tempDeme.simPheno=transect[i].getPheno(gps);
        tempDeme.ktable=transect[i].getK();
        simResults[i]=tempDeme;
    };
    return simResults;
};
