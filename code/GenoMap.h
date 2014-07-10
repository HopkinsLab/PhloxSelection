// This header file hold custom variable types that we use: including the genotype map (assuming two loci and independent assortment), haplotype and genotypes, and parent pairs, 

//  GenoMap.h
//  phloxsim
//
//  Created by Rafael Guerrero on 11/2/12.
//  Copyright (c) 2012 Robin Hopkins. All rights reserved.
//

#ifndef __phloxsim__GenoMap__
#define __phloxsim__GenoMap__

// General STL includes
#include <iostream>
#include <map>
using std::map;
#include <vector>
using std::vector;
#define NPHENOS 4

// Two inumerators for intensity and hue (with capital = 1 and lowercase =0)
enum intensity {
    y = 0,
    Y = 1
};

enum hue {
    h = 0,
    H = 1
};

// Haplotype holds the state of the loci (intensity and hue)
struct Haplotype{
    // Ideally, this would be the only struct to change if we want to add more loci
    // Apart from adding a variable to the struct, changes to the operators are necessary.
    
    intensity In;
    hue Hu;
    
    bool operator<(const Haplotype& other) const    {
		if(In <other.In) return true;
		if(In == other.In && Hu < other.Hu) return true;
		return false; }
    
	bool operator==(const Haplotype& other) const    {
		if(In == other.In && Hu == other.Hu) return true;
        return false; }

};

    
    
    
struct Genotype {
        // I would want the Genotype struct to be independent of the # of loci, so I added the Haplotype struct
        // HOWEVER For easy usage of the struct, two of the constructors are specific to 2 loci
        
        Haplotype M;
        Haplotype F;
        int phenotype;
        
        // Phenotype is 0-3 with phenotypes numbered in alphabetical order DB, DR, LR, LB
        
        Genotype(){}                    
        Genotype(const Haplotype D, const Haplotype E) {M = D; F = E; buildPhenotype();}
            // Constructor that takes two haplotypes (flexable to number of loci)
        Genotype(intensity Mi, hue Mh, intensity Fi, hue Fh) {
            // Constructor that takes two intensity and two hue types (h,H,Y,y)
            // Note that this constructor will need to be updated if a locus is added to the model
            Haplotype male = {Mi, Mh};
            Haplotype female = {Fi, Fh};
            M = male;
            F = female;
            buildPhenotype();
        }
        Genotype(const int Mi, const int Mh, const int Fi,const int Fh) {
            // Constructor that takes 0 & 1 as alleles (0,1,0,1)
            // Note that this constructor will need to be updated if a locus is added to the model
            Haplotype male = {(intensity)Mi, (hue) Mh};
            Haplotype female = {(intensity)Fi, (hue) Fh};
            M = male;
            F = female;
            buildPhenotype();
        }
    Genotype(const Genotype& other) {M=other.M; F=other.F; phenotype=other.phenotype;}
    
    bool recessiveIntensity(){
        // Test if homozygous recessive genotype (used for determining phenotype)
        if(M.In == F.In && M.In == y) return true;
        return false;
    };
   
    bool recessiveHue(){
        if(M.Hu == F.Hu && M.Hu == h) return true;
        return false;
    };
    void buildPhenotype(){
        // Assigns phenotype to Genotype struct
        if     (!recessiveIntensity() && !recessiveHue()) phenotype=0;      // Dark Blue phenotype
        else if(!recessiveIntensity() &&  recessiveHue()) phenotype=1;      // Dark Red phenotype
        else if( recessiveIntensity() && !recessiveHue()) phenotype=2;      // Light Blue phenotype
        else if( recessiveIntensity() &&  recessiveHue()) phenotype=3;      // Light Red phenotpye
        else{std::cout<<"Error in buildPhenotype:: undefined phenotype\n";}
    };
        bool operator<(const Genotype& other) const    {
            if(M <other.M) return true;
            if(M == other.M && F < other.F) return true;
            return false; }
        
        bool operator==(const Genotype& other) const    {
            if(M == other.M && F == other.F) return true;
            return false; }
    
};

    
    
typedef map<Genotype, int> gps_type;
typedef gps_type::iterator gps_i;
// Defines the parent pair struct 
    struct ParentPair{
        int m;                                         // Index of the male genotype in the frequency vector
        int f;                                         // Index of the female genotype in the frequency vector
        double fraction;                               // Fraction of offspring from that pair that produce a certain Genotype
        bool identity;                                 // Are the male and female phenotypes the same?
        int fePheno;                                   // Female phenotype in parent pair
        
        ParentPair(){}
        // Makes a parent pair using indexes from genotypes
        ParentPair( gps_i male,  gps_i fmale,  gps_i o){
            
            m=male->second;
            f=fmale->second;
            buildfraction(male->first, fmale->first, o->first);
            identity= male->first.phenotype == fmale->first.phenotype;
            fePheno= fmale->first.phenotype;
        }
    // Calculates the fraction of offspring genotype (o) from mating between male (m) and female (f) 
        void buildfraction(Genotype m, Genotype f, Genotype o){
            double maleProbInt=0;
            double maleProbHue=0;
            double fmaleProbInt=0;
            double fmaleProbHue=0;
            
            if(o.M.In == m.M.In) maleProbInt +=0.5; // if the offspring state is equal to one of the male's alleles, then
            if(o.M.In == m.F.In) maleProbInt +=0.5; // the prob of inheritance increases by 0.5
          
            if(o.M.Hu == m.M.Hu) maleProbHue +=0.5;
            if(o.M.Hu == m.F.Hu) maleProbHue +=0.5;
            
            if(o.F.In == f.M.In) fmaleProbInt +=0.5;
            if(o.F.In == f.F.In) fmaleProbInt +=0.5;
            
            if(o.F.Hu == f.M.Hu) fmaleProbHue +=0.5;
            if(o.F.Hu == f.F.Hu) fmaleProbHue +=0.5;
            
            fraction= maleProbInt*maleProbHue*fmaleProbInt*fmaleProbHue;
        }
    };
    
// Create phenoId. This is a list of genotypes that have each of the four phenotypes 
    struct PhenoID{
       vector<vector<int> > idList;
        PhenoID(){}
        PhenoID(gps_type& gps){
            idList.resize(NPHENOS);
            for (gps_i i=gps.begin(); i!=gps.end();++i){
                idList.at(i->first.phenotype).push_back(i->second);
                //std::cout<< "Genotype "<<i->second<< " has Phenotype "<<i->first.phenotype<<'\n';
            }
        
        }
    };
#endif /* defined(__phloxsim__GenoMap__) */


    
