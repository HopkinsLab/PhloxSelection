//
//  World.cpp
//  phloxsim
//
//  Created by Robin Hopkins on 10/23/12.
//  Copyright (c) 2012 Robin Hopkins. All rights reserved.
//
//#define FILESDBGR
//#define MIGDBGR
//#define PARAMDBGR
#include "World.h"

#include <iomanip>
#include <iostream>
using std::cout;
#include <vector>
using std::vector;
#include <string>
using std::string;
#include <fstream>
using std::ifstream;
#include <sstream>
using std::stringstream;
using std::istringstream;
#include<math.h>



World::World(Params* par, int key){
    // This function creates the world
    // It will initialize the genotype map (gps), 
    // the vector of demes (transect),
    // the vector of ParentPairs (parents),
    // and the neighborhood
    // For the construction of each deme we need: fitness coef (a pointer to a vector size 16), initial Genotype freqs, and an index number (nDeme)
    // That means this constructor will have to translate the coefs as input in the parameters (three per zone) and define which 
    // demes are in each zone
    a=par->a;
    makeGPS();
    build_wX(par);
    initialFreqs(readinput(par, key), par->totDemes);
    buildParents();
    buildTransect(par,key);
    buildNeighborhood(par->m, par->deltaDeme, par->migDis);
    
    
};

void World::makeGPS(){
    // Make genotype map, **note that it assumes two loci with two alleles each** (Mint, Mhue, Fint, Fhue)
    int counter=0;
    for(int i=0; i<2;++i){                              // i is male intensity
        for(int j=0; j<2;++j){                          // j is male hue
            for(int k=0; k<2;++k){                      // k is female intensity
                for(int l=0; l<2;++l){                  // l is female hue
                    Genotype dummyGeno(i,j,k,l);
                    gps[dummyGeno]=counter;
                    counter++;
                }
            }
        }
    }
    
    /*
    //Output map *** FOR DEBUGGING ONLY
  
     for (gps_i i=gps.begin(); i!=gps.end(); ++i) {
     Genotype dummyGeno = (*i).first;
     int dummyIndex = (*i).second;
     std::cout<<"Genotype: "<<dummyGeno.F.Hu<<dummyGeno.F.In<<dummyGeno.
         M.Hu<<dummyGeno.M.In<<" has index = "<<dummyIndex<<'\n';
     }
        */
    
};

void World::build_wX(Params* par){
    // Build the fitness vectors from the parameter file
    wX.resize(par->nZones);
    
    for (gps_i i=gps.begin(); i!=gps.end(); ++i) {
        Genotype dummyGeno = (*i).first;
        double fitness;
        
        for(int Z=0; Z<par->nZones;++Z){
            // For each of Z zones define the 16 fitness coeffecients for each of the genotypes based on the phenotype
            if      (!dummyGeno.recessiveHue() && !dummyGeno.recessiveIntensity()) fitness= par->w.at(Z)[0];
            else if ( dummyGeno.recessiveHue() &&  dummyGeno.recessiveIntensity()) fitness= par->w.at(Z)[2];
            else if ( dummyGeno.recessiveHue() && !dummyGeno.recessiveIntensity()) fitness= par->w.at(Z)[1];
            else if (!dummyGeno.recessiveHue() &&  dummyGeno.recessiveIntensity()) fitness= 1.0;           
            else { std::cout<<" Check build_wX function:: undefined fitness coefficient\n";}
            
            wX.at(Z).push_back(fitness);
           
            //FOR DEBUGGING- OUTPUT
            // std::cout<<"z "<<Z<<" index "<<index<<" w "<<wX.at(Z).at(index)<<'\n';
        }
    }
  
};
void World::buildTransect(Params* par, int key){
    // Create an instance of the PhenoID matrix
    FILESDBG("From World::buildTransect\n");
    PhenoID phenoGeno(gps);
    
    for(int i=0; i< par->totDemes; ++i){
        int Z=0;
        double borderNum= par->border[key] * par->totDemes;
        if( borderNum <= i) Z=1;
        vector<double>* wZone = & wX.at(Z);
        pX=startP[i];// defines initial genotype frequency for each deme
        FILESDBG("deme "<<i<<" genos: ");
                 for(int j=0; j<pX.size();++j){
                     FILESDBG(pX[j]<<" ");
                 }
        FILESDBG('\n');
        Deme d(i, wZone, pX, parents,phenoGeno);
        transect.push_back(d);
    };
    PARDBG("selection: ");
    for (int i=0;i<par->w.size();++i){
        for(int j=0; j<par->w[i].size(); ++j){
            double temp=par->w[i][j];
            PARDBG(temp<<" ");
             }
    }
    double tempm= par->m;
    double tempb= par->border[key];
    double tempa= par->a;
    PARDBG("m "<<tempm<<" bordernum "<<tempb<<'\n');
};

std::ostream& operator<<( std::ostream& os, const Genotype& gen )
{
    // this overloaded operator will make our output life easier, especially in debugging
    os<<gen.M.In <<gen.M.Hu <<gen.F.In <<gen.F.Hu;
    return os;
};

// Builds parent combintations "ParentPair" with males, females, and fraction of offspring
void World::buildParents(){
    parents.resize(gps.size());
    for (gps_i i=gps.begin(); i!=gps.end(); ++i) {                  // Goes through all male genotypes
        
        for (gps_i j=gps.begin(); j!=gps.end(); ++j) {              // Goes through all female genotypes
            
            for (gps_i k=gps.begin(); k!=gps.end(); ++k) {          // Goes through all offspring genotypes
                int index=(*k).second;
                ParentPair dummyP (i,j,k);
                if(dummyP.fraction !=0) {                           // Fraction is calculated inside construction of parent pair 
                    parents.at(index).push_back(dummyP);
                    
                // Output for checking
                   //std::cout<<"male "<<(*i).first<<", female "<< (*j).first<<", Offspring "<<(*k).first<<" index "<<(*k).second<<", fraction "<<dummyP.fraction<<" identity "<<dummyP.identity<<" male pheno "<<dummyP.phenoM<<'\n';
                     
                }
                
            }
        }
        
    }
    
};

// Build transect of Demes and declares the selection zones (assumed 2)

// For each deme identifies neighbors and the proportion of migrants contributed. We assume stepping stone model with migration from only imediate neighbors.
void World::buildNeighborhood(double& m, double& delta, int&migDis){
    vector<double>migvec;
    
    // calculate the migration probability from neighboring demes for the given variance (m) and space between demes (delta) with a ddesignated likelihood distribution (1=normal, 2=geometric)
    if (migDis==1){migvec=nbinMigRate(m, delta);}
    else if(migDis==2){migvec=gBinMigRate(m,delta);}
    else if(migDis==3){migvec=tenNormMigRate(m, delta);}
    else if(migDis==4){migvec=tenGBinMigRate(m,delta);}
    else if(migDis==5){migvec=tenMigraine(m);}
    else {migvec=nbinMigRate(m,delta);
        cerr<<"migration distribution might be wrong. check";
    };
    vector < vector <double> > mig_prob;
    mig_prob.clear();
	int nDemes=transect.size();
    mig_prob.resize(nDemes);
	
	for(int i = 0; i < nDemes; i++) mig_prob.at(i).resize(nDemes); // make neighbor vector for each deme size of number of demes
	
    for(int i = 0; i < nDemes; i++){        // i is the focal deme - go through each focal deme and add neighbors
        double migtot=0;
        
        for(int mbin = 0; mbin < migvec.size(); mbin++){
            MIGDBG("bin "<<mbin<<" ");
            //go through the bins in the migration probability vector and assign RIGHT neighbors migration rates
            int nbRT=i+mbin;// the right side neighbor that lines up to the mbin in the migration vector
            
            if(nbRT<nDemes){
                mig_prob.at(i).at(nbRT)+= migvec[mbin];		// proportion of migrats from neighbor (nb) that is mbin demes away from focal deme (i)
                migtot+=migvec[mbin];
                MIGDBG("RTmigamount "<<mig_prob.at(i).at(nbRT)<<" for deme"<<i<<" goes to deme "<<nbRT<<" ");
            }
            // if nb (the needed neighbor) is beyond the number of demes, this puts the mig prob in the correct deme
            else {
                MIGDBG("else ");
                double bin=nDemes-(nbRT-nDemes+1);
                
                mig_prob.at(i).at(bin)+=migvec[mbin];
                migtot+=migvec[mbin];
                MIGDBG("RTmigamount "<<mig_prob.at(i).at(bin)<<" for deme"<<i<<" goes to deme "<<bin<<" ");
            }
            //go through the bins in the migration probability vector and assign LEFT neighbors migration rates
            int nbLT=i-mbin;            //the left side neighbor that lines up to the mbin in the migration vector
            if(nbLT>=0){
                mig_prob.at(i).at(nbLT)+= migvec[mbin];		// proportion of migrats from neighbor (nb) that is mbin demes away from focal deme (i)
                migtot+=migvec[mbin];
                MIGDBG("LTmigamount "<<mig_prob.at(i).at(nbLT)<<" for deme"<<i<<" goes to deme "<<nbLT<<'\n');
            }
            //if nbLT is to the left of deme zero this puts the migration in the correct deme
            else {
                double binlft=-(nbLT)-1;
                mig_prob.at(i).at(binlft)+=migvec[mbin];
                migtot+=migvec[mbin];
                MIGDBG("LTmigamount "<<mig_prob.at(i).at(binlft)<<" for deme"<<i<<" goes to deme "<<binlft<<'\n');
            }
        };
        
        if(migtot<1-ERR||migtot>1+ERR){
            std::cerr<<"Migration rate is wrong. exiting "<<migtot;
            exit(1);
        }
        
        MIGDBG("migration total for deme"<<i<<": "<<migtot<<'\n');
    };
    // Fill in the Neighbors (with non-zero migration) for each Deme in the transect
    for(int i = 0; i < nDemes; i++){
        MIGDBG("Deme "<<i<<" has migration ");
        for(int j = 0; j < nDemes; j++){
            if (mig_prob.at(i).at(j)!=0) {
                Neighbor d(&transect.at(j), mig_prob.at(i).at(j));     //Creates a neighbor d (pointer to deme giving migrants and fraction of migrants)
                transect.at(i).putNeighbor(d);                          // Puts neighbor d into a list of neighbors for Deme i.
                
                MIGDBG("from deme "<<j<<" "<<mig_prob.at(i).at(j)<<" ");
            }
        }
        MIGDBG('\n');
    }
    
};




