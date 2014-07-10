//  The header file for the Deme class. Contains the core functions that happen within (and between) demes - selection, migration, mating.
//  Deme.h
//  phloxsim
//
//  Robin Hopkins on 11/2/12.

//

#ifndef __phloxsim__Deme__
#define __phloxsim__Deme__
#define ER 0.001
// Includes from STL
#include <iostream>
#include <vector>
using std::vector;

// This is an include from Boost library: boost.org 
#include <boost/shared_ptr.hpp>
using boost::shared_ptr;

// Our program specific includes
#include "GenoMap.h"

// Forward decleration of types used by Demes.
class World;
struct Neighbor;
struct ParentPtr;


class Deme{
   
public:
    Deme(unsigned int nDemeX, vector<double>* wATzone, vector<double> pX, const vector<vector<ParentPair> >& parents, PhenoID& phenoGeno);     // Deme constructor
    ~Deme();    // Deme destructor
    
    //functions:
    void buildParentPtrs();
    void selectD();                                             // Declare function selection
    void migrateD();                                            // Declare migration function acting on males
    void matingD(vector<vector< ParentPair> >& parents, gps_type& gps, double a);        // Declare function mating
    void makeKD(vector<vector< ParentPair> >& parents, gps_type& gps, double a);// Make vector of vectors
    
    
    vector<double> getPheno(gps_type& gps);
    vector<shared_ptr<double> > getP();                                  
    double getP(Genotype geno);
    double getP(int geno);
    
    vector< shared_ptr<double> > getPS();
    
    double getPS(int geno);
    
    vector<shared_ptr<double> > getPSptr();                            // Pointer to frequencies after selection
    
    
    vector<shared_ptr<double> > getPM();
    double getPM(int geno);
    
    double getPDel();
    
    vector<vector<double> > getK();
    
    void putP(vector<shared_ptr<double> > pX);
    void putW(vector<double> * wATzone);
    void putNDeme(unsigned int nDemeX);
    void putNeighbor(Neighbor newneighbor);
    double assort(gps_type& gps);
    
private:
    // Data, need to be given at construction or with buildNeighborhood
    unsigned int nDeme;                                         // Location in transect of deme
    vector<shared_ptr<double> > p;                              // Starting genotype frequencies
    vector<double>* w;                                          // Pointer to selection coefficients
    vector<Neighbor>n;                                          // Neighbors for migration
    vector<vector< ParentPtr> > fams;                           // Vector of vectors of parent pair structs
    
    // Calculated inside of Deme
    vector< vector< shared_ptr<double> > > phenoGenoPtr;
    vector<shared_ptr<double> > pS;                                //  allele frequencies after selection
    vector<shared_ptr<double> > pM;                             // Male allele frequencies after migration
    double pDel;                                                // Max dif across genos in freq from start to end of gen
    vector<vector<double> > k;                                  // Probability of offspring from female phenotype(with first index ) having phenotype (second index);

};

struct Neighbor{
    // List of pointers to neighbors giving migrants and the migration rates.
    Deme* nPoint;
    double m;
    
    Neighbor(){}
    ~Neighbor(){}
    Neighbor(Deme* x,double mX):nPoint(x),m(mX){}
};
// Beware - this is NOT the parent pair.  This takes the parent pair indexes (from ParentPair) and makes pointers to the genotype frequencies specific to each deme. This actually points to the p (females) and pM (Males)
struct ParentPtr{
    shared_ptr<double> mfreq;               // Pointer to the male genotype frequency in pM vector
    shared_ptr<double> ffreq;               // Pointer to the female genotype frequency in p vector
    double fraction;                        // Fraction of offspring from that pair that produce a certain Genotype
    bool identity;                          // Are the male and female phenotypes the same?
    int fPheno;                             // female phenotype
    
    
    ParentPtr();
    ParentPtr( shared_ptr<double> mo, shared_ptr<double> fo, double frac, bool i, int fP){
        mfreq = mo;
        ffreq = fo;
        fraction= frac;
        identity= i;
        fPheno =fP;
        }   
};


  
#endif /* defined(__phloxsim__Deme__) */

