//
//  utilities.cpp
//  phloxsim
//
//  Created by Robin Hopkins on 3/12/13.
//  Copyright (c) 2013 Robin Hopkins. All rights reserved.
//
//#define MIGDBGR
#include "utilities.h"
#define RES 0.001    //amount of residual area remaining in migration distribution area 

//Interpolate function
double interpolate (double left, double right, double gain, double x_axis){
    double slope = (right -left)/ x_axis;
    double val = slope *gain;
    return left + val;
};

// PDF for standard normal distribution.
double pdf (double dis){
    double out;
    out=1/(sqrt(2*M_PI))*exp(-(.5)*dis*dis);
    return out;
};
// aproximation of CDF for normal distribution Zelen & Severo (1964)
double normCDF (double dis, double Mvar){
    double disSD=dis/(sqrt(Mvar));
    const double b0 = 0.2316419, b1 = 0.319381530, b2 = (0-0.356563782), b3 = 1.781477937, b4 = (0-1.821255978), b5 = 1.330274429;
    double t=1/(1+b0*disSD);
    double t2=pow(t,2), t3=pow(t,3), t4=pow(t,4), t5=pow(t,5);
    double temppdf=pdf(disSD);
    double cdf=1-temppdf*(b1*t + b2*t2 + b3*t3 + b4*t4 + b5*t5);
    //UTDBG("cdf for dis "<<dis<<" is "<<cdf<<" ");
    return cdf;
};
// creates a vector with migration probabilities for bins starting with not migrating [0] to right neighbor [1]
// We are dividing our continuous dispersal function into "bins" that correspond to our demes.  There is migration between these bins/demes proportional to the area under the dispersal curve in the bin.  So first we determin our bins (left side and right side of deme) and then calculate the area under the dispersal curve in that bin.  we only do this with enough bins to get most of the distribution (less ER).  We then relativize the bins so they add to one. 
  // takes sigma squared (Mvar) and the deme delta - distance between demes in km) gives a vector of migration probabilities given that sigma squared is in kmsqured.


vector <double> nbinMigRate(double Mvar, double delta){
    MIGDBG("Mvar "<<Mvar<<" and deltademe "<<delta<<"\n");
    vector<double> migRateVec;
    double migtot=0.0;
    int bin=0;
    
    while (migtot<1-RES){
        double binright = (bin+0.5)*delta;              // this is the right border of the bin (1/2 way to the neighbor deme
        double binleft = (bin-0.5)*delta;               // this is the left border of the bin (1/2 way to the next deme)
        double tempmig=0;
        double templeft=0;
        double tempright=normCDF(binright, Mvar);                // this is the cdf from the right border of the bin
        
        if (bin==0){templeft=0.5;}
        else {templeft=normCDF(binleft, Mvar);}          //this is the cdf from the left border of the bin          
                
        tempmig=tempright-templeft; // this is the area holding the probability of migrating to the bin (area between left and right border)
        migtot+=2*tempmig;          // migration is symmetric so the area counts twice toward the total probabilty of migration
        migRateVec.push_back(tempmig);
        MIGDBG("migbin "<<bin<<" "<<tempmig<<" migtot "<<migtot<<" \n");
        
        ++bin;
    };
    double migtotal=0.0;
        MIGDBG("final mig rate for ");
    for (int i=0;i<migRateVec.size();++i){
        migRateVec[i]=migRateVec[i]/migtot;
        migtotal+=migRateVec[i];
        MIGDBG("bin"<<i<<": "<<migRateVec[i]<<" ");
    };
        MIGDBG('\n');
    if(migtotal<0.5-ERR||migtotal>0.5+ERR){std::cerr<<"check binMigRate. total migration is"<<migtotal<<"\n";exit(1);}
    
    return migRateVec;
};



vector <double> tenNormMigRate(double Mvar, double delta){
    MIGDBG("Mvar "<<Mvar<<" and deltademe "<<delta<<"\n");
    vector<double> migRateVec;
    double migtot=0.0;
        
    for(int bin=0; bin<10;++bin){
        double binright = (bin+0.5)*delta;              // this is the right border of the bin (1/2 way to the neighbor deme
        double binleft = (bin-0.5)*delta;               // this is the left border of the bin (1/2 way to the next deme)
        double tempmig=0;
        double templeft=0;
        double tempright=normCDF(binright, Mvar);                // this is the cdf from the right border of the bin
        
        if (bin==0){templeft=0.5;}
        else {templeft=normCDF(binleft, Mvar);}          //this is the cdf from the left border of the bin
        
        tempmig=tempright-templeft; // this is the area holding the probability of migrating to the bin (area between left and right border)
        migtot+=2*tempmig;          // migration is symmetric so the area counts twice toward the total probabilty of migration
        migRateVec.push_back(tempmig);
        MIGDBG("migbin "<<bin<<" "<<tempmig<<" migtot "<<migtot<<" \n");
       
    };
    double migtotal=0.0;
    MIGDBG("final mig rate for ");
    for (int i=0;i<migRateVec.size();++i){
        migRateVec[i]=migRateVec[i]/migtot;
        migtotal+=migRateVec[i];
        MIGDBG("bin"<<i<<": "<<migRateVec[i]<<" ");
    };
    MIGDBG('\n');
    if(migtotal<0.5-ERR||migtotal>0.5+ERR){std::cerr<<"check binMigRate. total migration is"<<migtotal<<"\n";exit(1);}
    
    return migRateVec;
};




//PDF for geometric distribution
double geoPDF(double dis, double mVar, double deltaX){    
    
    double gpdf=0;          //probability of dispersal using the geometric distribution
    double dXsqr=deltaX*deltaX;
    double sqtemp=sqrt(dXsqr+8*mVar);  //just to get squareroot
    double p=0.0;
    //The mVar is the variance in the geometric distribution, with this varience, p can be calculated.
    if (dXsqr==mVar){p=2/3;}
    else {p=(3*dXsqr-deltaX*sqtemp)/(2*(dXsqr-mVar));}               
    
    double temp=pow(1-p,dis);
    if(dis==0){gpdf=p/2;}   // Divide by 2 so that the migration for loop works for left and right migration and so the making of the bins stops correctly.
    else {gpdf=p*temp/2;}
    
    return (gpdf);
}
// Creates a vector of migrations if migration likelihood is a geometric distribution. The same as above but for a different distribution.
vector<double> gBinMigRate(double Mvar, double delta){
    MIGDBG("Mvar "<<Mvar<<" and deltademe "<<delta<<"\n");
    vector<double> migRateVec;
    double migtot=0.0;
    int bin=0;
    
    while (migtot<1-RES){
        double tempmig=0;
        tempmig=geoPDF(bin,Mvar,delta);
        migtot+=2*tempmig;
              
        migRateVec.push_back(tempmig);
        MIGDBG("migbin "<<bin<<" "<<tempmig<<" migtot "<<migtot<<" \n");
        
        ++bin;
    };
    double migtotal=0.0;
    MIGDBG("final mig rate for ");
    for (int i=0;i<migRateVec.size();++i){
        migRateVec[i]=migRateVec[i]/migtot;
        migtotal+=migRateVec[i];
        
        MIGDBG("bin"<<i<<": "<<migRateVec[i]<<" ");
    };
    MIGDBG('\n');
    if(migtotal<0.5-ERR||migtotal>0.5+ERR){std::cerr<<"check binMigRate. total migration is"<<migtotal<<"\n";exit(1);}
    
    return migRateVec;
};

vector<double> tenGBinMigRate(double Mvar, double delta){
    MIGDBG("Mvar "<<Mvar<<" and deltademe "<<delta<<"\n");
    vector<double> migRateVec;
    double migtot=0.0;
    
    for(int bin=0; bin<10; ++bin){
        double tempmig=0;
        tempmig=geoPDF(bin,Mvar,delta);
        migtot+=2*tempmig;
        
        migRateVec.push_back(tempmig);
        MIGDBG("migbin "<<bin<<" "<<tempmig<<" migtot "<<migtot<<" \n");
        
    };
    double migtotal=0.0;
    MIGDBG("final mig rate for ");
    for (int i=0;i<migRateVec.size();++i){
        migRateVec[i]=migRateVec[i]/migtot;
        migtotal+=migRateVec[i];
        
        MIGDBG("bin"<<i<<": "<<migRateVec[i]<<" ");
    };
    MIGDBG('\n');
    if(migtotal<0.5-ERR||migtotal>0.5+ERR){std::cerr<<"check binMigRate. total migration is"<<migtotal<<"\n";exit(1);}
    
    return migRateVec;
};

vector<double> tenMigraine(double Mvar){
    vector<double>migRateVec;
    double migtot=0.0;
    
    for(int bin=0; bin<10; ++bin){
        double tempmig;
        if(bin==0){ tempmig=(1-Mvar)/2;}
        else{ tempmig=(Mvar/18);}
        migtot+=2*tempmig;
        migRateVec.push_back(tempmig);
    }
    
    if(migtot<1-ERR||migtot>1+ERR){std::cerr<<"check tenmigraine bin rates total migration is "<<migtot<<"\n";exit(1);}
    return migRateVec;
};


//Prior function. determines the prior distribution

double calcPrior(double param,double priorRate, double priorMax){
    double s=priorMax;
    double r=priorRate;
    double priorP=0.0;
    double temppri=exp(r*(param-s))/(exp(r*(param-s))+1);
    priorP=log(temppri);
    return priorP;
};


//Output to file function
void outToFile (vector<DemeOut> &s, const char* ss){
    
   
ofstream fout( ss,ios::out);  // output file stream
	fout.precision(8);
    
    fout<<"Deme dbFreq drFreq lbFreq lrFreq DBdb DBdr DBlb DBlr DRdb DRdr DRlb DRlr LBdb LBdr LBlb LBlr LRdb LRdr LRlb LRlr \n";
    for(int i=0; i<s.size();++i){
        fout<<i<<" ";
        for(int j=0; j<s[i].simPheno.size();++j){
            fout<<s[i].simPheno[j]<<" ";
        };
        for(int k=0; k<s[i].ktable.size(); ++k){
            for (int j=0; j<s[i].ktable[k].size(); ++j){
                fout<<s[i].ktable[k][j]<<" ";
            };
        };
        fout<<'\n';
    };
    
    fout.close();
};
vector<double> NR_to_vec(VecDoub_I inVec){
    vector<double> outVec;
    for(int i=0;i<inVec.size();++i) outVec.push_back(inVec[i]);
    return outVec;
};

VecDoub vec_to_NR(vector<double> inVec){
    VecDoub outVec;
    outVec.resize(inVec.size());
    for(int i=0;i<inVec.size();++i) outVec[i]=inVec[i];
    return outVec;
};

double inverseLogit(double a, vector<double> range){
    //take any number a (-inf,inf) and return a value within range
    double p = exp(a)/(1.0+exp(a));//p is a number in [0,1)
    double length= (range[1]-range[0]);//Scaling to the par of interest
    double out=p*length + range[0];
    
    return out;
};

double Logit(double p, vector<double> range){
    //take any p within range and return a (-inf,inf)
    double ptemp= (p-range[0])/(range[1]-range[0]); //setting ptemp~[0,1]
    double a = log(ptemp)-log(1-ptemp); // this is the actual transform, a is (-inf,inf)
    if(ptemp<0){std::cerr<<"Error in Logit():: cannot do Log of negatives\n";exit(1);}
    return a;
};

std::ostream& operator<<(std::ostream& os, const vector<double>& g){ for (int i=0;i < g.size();++i) std::cout<<g.at(i)<<" "; return os;};

std::ofstream& operator<<(std::ofstream& os, const vector<double>& g){ for (int i=0;i < g.size();++i) std::cout<<g.at(i)<<" "; return os;};
