//
//  utilities.h
//  phloxsim
//
//  Created by Robin Hopkins on 3/12/13.
//  Copyright (c) 2013 Robin Hopkins. All rights reserved.
//

#ifndef __phloxsim__utilities__
#define __phloxsim__utilities__

#include <vector>
using std::vector;

#include <iostream>
using std::ios;
using std::cout;
#define _USE_MATH_DEFINES
#include<math.h>
#include "nr3.h"

#include <fstream>
using std::ifstream;
using std::ofstream;

#include <sstream>
using std::stringstream;
#define ERR 0.0001


//Struct that describes the simulated output from each deme. Includes phenotype frequencies and ktables
struct DemeOut{
    vector<double> simPheno;                            // Expected phenotype frequencies
    vector<vector<double> > ktable;                     // Expected phenotype of offspring given mom phenotype
    
    DemeOut(vector<double> sP,vector<vector<double> > kt):simPheno(sP),ktable(kt){}
    DemeOut(){};
    
};

// Phenotype data from field collected mom's and offspring
struct FamData{
    vector<int> offCount;                       // Counts of offspring phenotype for a mom
    int momPheno;                               // Mom phenotype (0-3)
};

// Phenotype counts from field populations
struct ObsData{
   // double rawPos;                             // position from start of transect
    double position;                          // Position from start of transect
    vector<int> PhenoCounts;                  // Counts of phenotypes from a populations
    vector<FamData> momOff;                   // The vector of family data from that population
};


//Interpolate function
double interpolate (double left, double right, double gain, double x_axis);

// aproximation of CDF for normal distribution
double normCDF (double dis, double Mvar);
// creates a vector with probabilities of migrating using NORMAL distribution
vector <double> nbinMigRate(double Mvar, double delta);
vector <double> tenNormMigRate(double Mvar, double delta);

//CDF for geometric distribution
double geoCDF(double dis, double mVar, double deltaX);
// creates a vector with probabilities of migrating using a geometric distribution

vector<double> gBinMigRate(double Mvar, double delta);
vector <double> tenGBinMigRate(double Mvar, double delta);
vector<double> tenMigraine(double Mvar);


double calcPrior(double param,double priorRate, double priorMax);


//Output to file function
void outToFile (vector<DemeOut> &s, const char* ss);

vector<double> NR_to_vec(VecDoub_I inVec);

VecDoub vec_to_NR(vector<double> inVec);

double inverseLogit(double a, vector<double> range);

double Logit(double p, vector<double> range);

std::ostream& operator<<(std::ostream& os, const vector<double>& g);
std::ofstream& operator<<(std::ofstream& os, const vector<double>& g);

#endif /* defined(__phloxsim__utilities__) */
#ifdef MIGDBGR
#define MIGDBG(x) std::cout << x
#else
#define MIGDBG(x)
#endif

