// The deme contains all the meat of the simulation. Defines all the functions (selection, migration, mating)

//  Deme.cpp
//  phloxsim
//
//  Created by Robin Hopkins on 11/2/12.
//  Copyright (c) 2012 Robin Hopkins. All rights reserved.
//


#include "Deme.h"


#include <iostream>
#include <math.h>
#include <boost/foreach.hpp>

#include "Params.h"




Deme::Deme(unsigned int nDemeX, vector<double>* wATzone, vector<double > pX, const vector<vector<ParentPair> >& parents, PhenoID& phenoGeno){                                // Deme constructor
    nDeme = nDemeX;                                     // Number of demes
    w = wATzone;                                        // Pointer to fitness zone
    //pS.resize(pX.size());                               // Selection vector
    
    for(int i=0;i<pX.size();++i){                       // Initializes vector of pointers for gene freq (p and pM)
        double ddd=pX.at(i);
        shared_ptr<double> forP(new double(ddd));
        p.push_back(forP);
        shared_ptr<double> forPs(new double(0.0));
        pS.push_back(forPs);
        shared_ptr<double> forPM(new double(0.0));
        pM.push_back(forPM);
    }
    
  
    for(int i=0; i<parents.size();++i){                 // Make vector of to parentPtr's to create offspring 
        fams.resize(parents.size());
        for(int j=0; j<parents.at(i).size();++j){
            int dad =parents.at(i).at(j).m;             // Takes the index from the ParentPair vector (world) and uses it to point at the correct male genotype frequency.
            int mom = parents.at(i).at(j).f;
            double fraction = parents.at(i).at(j).fraction;
            bool id = parents.at(i).at(j).identity;
            int fePhe=parents.at(i).at(j).fePheno;;
          ParentPtr dummy(pM.at(dad), pS.at(mom), fraction, id,fePhe);
            fams.at(i).push_back(dummy);
            
            // OUTPUT for debugging lists all the female frequencies in all the parentPtrs
            //std::cout<<"mom "<<mom<< " fepheno "<<fePhe<<'\n';
        }
    }
    
    // This is a vector of vectors of pointers to male genotype frequencies after migration. It is created from the genotypes listed in the phenoGeno struct (see GenoMap).  This is used to calculate the total assortative mating freq.
    phenoGenoPtr.resize(NPHENOS);
    for(int i=0; i<phenoGenoPtr.size();++i){
        for(int j=0; j<phenoGeno.idList.at(i).size();++j){
            phenoGenoPtr.at(i).push_back(pM.at(phenoGeno.idList.at(i).at(j)));
            //std::cout<< "phenotpe "<<i<< " with genotype "<<j<<'\n';
        }
    }
    
    };

Deme::~Deme(){                                                  // Deme destructor
    
};


vector<shared_ptr<double> > Deme::getP(){                                    // Get vector of genotype frequencies
    return p;
};

double Deme::getP(int geno){                                    // Get frequency for a given genotype
    return *p[geno];
};


vector< shared_ptr<double> > Deme::getPS(){                                    // Get vector of male genotype frequencies after selection
    return pS;
};

double  Deme::getPS(int geno){                                    // Get frequency for a given male genotype after selection
    return *pS[geno];
};

vector<shared_ptr<double> > Deme::getPM(){                       // Get vector of male genotype frequencies after migration
    return pM;
};

double Deme::getPM(int geno){                                    // Get frequency for a given male genotype after migration
    return *pM[geno];
};

double Deme::getPDel(){                              // Get max of genotype frequencies differences between two generations
    return pDel;
};

vector<vector<double> > Deme::getK(){                       // Get the K offspring probability table
    return k;
    
};

vector<shared_ptr<double> > Deme::getPSptr(){                         // Get vector of male genotype frequencies after selection
    return pS;
};



void Deme::buildParentPtrs(){
    
    
};



// Does selection on all genotype in this deme
void Deme::selectD(){
    double wB=0.0;                                              // Mean fitness
    for(int i=0; i<p.size(); ++i) {
        *pS[i]=(*p[i]) * w->at(i);
        wB+=*pS[i];
    };    
    for(int i=0; i<p.size(); ++i) {
        *pS[i]/=wB;                                              // Male genotype frequency after selection
        
    };
    
    //OUTPUT CHECK
    //std::cout<< "w bar is "<<wB<<'\n';
};

// Migration from neighbors into this deme
void Deme::migrateD(){
      for (int i=0; i<p.size(); ++i){       // Loop through each genotype
          double count=0.0;                 // New frequency of males in deme after migration
          double migtot=0;
          
          // Total mig rate of all other demes into this deme
          for(int j=0; j<n.size();++j){
              count +=(n.at(j).nPoint->getPS(i) * n.at(j).m);
              migtot+=n.at(j).m;
          }
           *pM[i]=count;
          if(migtot<1-ER||migtot>1+ER){
              std::cerr<<"Migration rate is wrong in deme. "<<migtot<<" exiting";
              exit(1);
          }
      }
};

//Total frequency of matings that are assortative by phenotype. goes through each genotype and multiplies the female frequency by each male frequency with the same phenotype (as listed in the phenoGenoPtr)
double Deme::assort(gps_type& gps){
    double atemp=0.0;
      
    for (gps_i j=gps.begin(); j!=gps.end();++j){
        int Fnum=(*j).second;
        Genotype Fgeno=(*j).first;
        int Fpheno=Fgeno.phenotype;
                
        for (int k=0; k<phenoGenoPtr.at(Fpheno).size();++k){
             atemp+=*p.at(Fnum) * *phenoGenoPtr.at(Fpheno).at(k);
        }
             }
    return atemp;
};




// Does recombination and mating and produced offspring genotype frequencies for next generation
void Deme::matingD(vector<vector< ParentPair> >& parents, gps_type& gps, double a){
    //a=0.0;
    double atot=assort(gps);
    //std::cout<< "total assortative matings "<<atot<<'\n';
    pDel=0;
    vector<double>tempfreq;
        
    for (gps_i j=gps.begin(); j!=gps.end(); ++j) {
        double counter=0.0;
        tempfreq.resize(fams.size());
        int index =(*j).second;
       // std::cout<<"Index "<<index<<", fams "<<fams.at(index).size()<<'\n';
    
        for (int i=0; i < fams[index].size(); ++i){
            
            //std::cout<<"index "<<index<<" pair # "<<i<<" male freq "<<(*fams[index][i].mfreq) <<" female freq " <<*fams[index][i].ffreq<<'\n';
            
            // Random mating that creates offspring genotype
            counter += (*fams[index][i].mfreq) *(1-a)* (*fams[index][i].ffreq) * fams[index][i].fraction;
           // std::cout<<"a "<<a<<'\n';
          //std::cout<<counter<<'\n';
            // Assortative mating that creates offspring genotype
           if(fams[index][i].identity)
               counter += a * (*fams[index][i].mfreq) * (*fams[index][i].ffreq) * fams[index][i].fraction / atot;
            
        }
        
        // Identifies max change in genotype frequency
        pDel=std::max(fabs(counter - (*p.at(index))), pDel);
        
        tempfreq.at(index)= counter;
        
        
    }
         //std::cout<<"a equals "<<a<<'\n';
    for (int i=0; i < tempfreq.size(); ++i){
        *p.at(i)=tempfreq[i];
    }
};
// Calculates probablity of a female of each phenotype having an offspring of each phenotype (in a 4X4) matrix

void Deme::makeKD(vector<vector< ParentPair> >& parents, gps_type& gps, double a){
    //a=0.0;
    double atot=assort(gps);
    vector<double> kTot;                                // Total fraction of offspring from a given female pheno
    kTot.resize(NPHENOS);
    
    k.clear(); k.resize(NPHENOS); //RG V-5-13: added clear() -- bug fix 
    
    double tot=0.0;
    for (int index=0; index!=k.size();++index){
        k[index].resize(NPHENOS);
    }
    
    for (gps_i j=gps.begin(); j!=gps.end(); ++j) {      //Loop through each genotype
        int index =(*j).second;
        int offpheno=(*j).first.phenotype;              // Phenotype of offspring
        
        for (int i=0; i < fams[index].size(); ++i){     // Loop through each possible mating pair
            double counter=0.0;
            int fpheno=fams[index][i].fPheno;           // Phenotype of female
            // Mating equation to calculate the fraction of offspring from each mating
            counter += (*fams[index][i].mfreq) *(1-a)* (*fams[index][i].ffreq) * fams[index][i].fraction;
            
            if(fams[index][i].identity)
                counter += a * (*fams[index][i].mfreq) * (*fams[index][i].ffreq) * fams[index][i].fraction / atot;
            
            k[fpheno][offpheno]+= counter;              // Add the fraction to the k matrix
            kTot[fpheno]+=counter;
            tot+=counter;
           // if (fpheno==3) std::cout<<"offspring "<<offpheno<<" tot prob "<<k[fpheno][offpheno]<<" counter "<<counter<<'\n';
        }
        
    }
   // std::cout<< "tot "<<tot<<" pheno0 "<<kTot[0]<<" pheno1 "<<kTot[1]<<" pheno2 "<<kTot[2]<<" pheno3 "<<kTot[3]<<'\n';
    
    for(int fe=0; fe<k.size(); ++fe){
       for(int of=0; of<k.at(fe).size();++of){
           if (kTot[fe]!=0)
           k[fe][of]/=kTot[fe];
        }
    }
};


vector<double> Deme::getPheno(gps_type& gps){
    vector<double> phenFreq;
    phenFreq.resize(NPHENOS);
    for(gps_i i=gps.begin(); i!=gps.end(); ++i){
        int pheno=(*i).first.phenotype;
        int index=(*i).second;
        phenFreq[pheno]+=*p[index];
    }
    return phenFreq;
};

void Deme::putP(vector<shared_ptr<double> > pX){                // Put initial frequencies in p
    p = pX;
};


void Deme::putW(vector<double> * wATzone){                      // Put selection coefficients in w
    w = wATzone;
};


void Deme::putNDeme(unsigned int nDemeX){                       // Assign deme number to deme
    nDeme = nDemeX;
};

void Deme::putNeighbor(Neighbor newneighbor){                   // Adds neighbors to n, the vector of neighbors
    n.push_back(newneighbor);
};
