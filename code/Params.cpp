//
//  paramS.cpp
//  phloxsim
//
//  Created by Robin Hopkins on 10/19/12.
//  Copyright (c) 2012 Robin Hopkins. All rights reserved.
//
#define PDBGR
#include "Params.h"

#include <string>
using std::string;
#include <fstream>
using std::ifstream;
#include <sstream>
using std::istringstream;
#include <vector>
using std::vector;
#include <math.h>


#include "utilities.h"
#include"nr3.h"

Params::Params(const char* inputfile){
    
    ifstream fin( inputfile);
    
	string line;
	int temp;
	double doubtemp;
    
    getline(fin, line);
    istringstream iss(line);
    
    
 	while(iss >> temp){
		totDemes=temp;
	}
    PDBG("Number of Demes = "<<totDemes)
    
    getline( fin, line);
	iss.clear();
	iss.str(line);
    while(iss >> temp){
		nZones=temp;
	}
    PDBG("Number of Zones = "<<nZones)
    
    border.clear();
    getline( fin, line);
        iss.clear();
        iss.str(line);
        while(iss >> doubtemp){
            border.push_back(doubtemp);
            //PDBG("world"<<i<<" border = "<<border[i])
    
    };
    
    
    
    getline( fin, line);
	iss.clear();
	iss.str(line);
    while(iss >> doubtemp){
		m=doubtemp;
	}
    PDBG("Migration prob = "<<m)
    
//    getline( fin, line);
//	iss.clear();
//	iss.str(line);
//    while(iss >> doubtemp){
//		a=doubtemp;
//	}
//    PDBG("Assortative mating = "<<a)
//    
    getline( fin, line);
	iss.clear();
	iss.str(line);
    while(iss >> temp){
		maxGens=temp;
	}
    PDBG("Max number of generations = "<<maxGens)
    
    getline( fin, line);
	iss.clear();
	iss.str(line);
    while(iss >> doubtemp){
		maxChange=doubtemp;
	}
    if(maxChange==0) maxChange= 0.000000000001;
    PDBG("Max frequency change = "<<maxChange)
    
    w.resize(nZones);
    for(int i=0; i<w.size(); ++i){
        
        getline( fin, line);
        iss.clear();
        iss.str(line);
        while(iss >> doubtemp){
            w.at(i).push_back(doubtemp);
        }
        PDBG("Relative fitness in zone "<<i<<": DB="<<w.at(i).at(0)<<", DR="<<w.at(i).at(1)<<", LR="<<w.at(i).at(2))
    }
    
};

Params::Params(){
    
    ifstream inNuissance("inPhlox_nuissance.params");
    if(inNuissance.is_open()){
        PDBG("Reading 'inPhlox_nuissance.params':")
        string line;
        int temp;
        double doubtemp;
        
        getline( inNuissance, line);
        istringstream iss(line);
        while(iss >> temp){
            nWorlds=temp;
        }
        PDBG("Number of worlds = " <<nWorlds)
        
        getline( inNuissance, line);
        iss.clear();
        iss.str(line);        
        while(iss >> doubtemp){
            wName.push_back(doubtemp);
        }
        PDBG("Name of first world "<<wName[0]);
        if (nWorlds != wName.size()){std::cerr<<"Error in inPhlox_nuissance check number and names of worlds\n";exit(1);}
        getline( inNuissance, line);
        iss.clear();
        iss.str(line);
        
        while(iss >> temp){
            totDemes=temp;
        }
        PDBG("Number of Demes = "<<totDemes)
        
        getline( inNuissance, line);
        iss.clear();
        iss.str(line);
        while(iss >> temp){
            nZones=temp;
        }
        PDBG("Number of Zones = "<<nZones)
        
        getline( inNuissance, line);
        iss.clear();
        iss.str(line);
        while(iss >> temp){
            maxGens=temp;
        }
        PDBG("Max number of generations = "<<maxGens)
        
        getline( inNuissance, line);
        iss.clear();
        iss.str(line);
        while(iss >> doubtemp){
            maxChange=doubtemp;
        }
        if(maxChange==0) maxChange= 0.000000000001;
        PDBG("Max frequency change = "<<maxChange)
        
        getline( inNuissance, line);
        iss.clear();
        iss.str(line);
        while(iss >> doubtemp){
            temperature=doubtemp;
        }
        PDBG("Starting Temperature = "<<temperature)
        
        getline( inNuissance, line);
        iss.clear();
        iss.str(line);
        while(iss >> doubtemp){
            temp_step=doubtemp;
        }
        PDBG("Temperature step= "<<temp_step)

        getline( inNuissance, line);
        iss.clear();
        iss.str(line);
        while(iss >> doubtemp){
            tolerance=doubtemp;
        }
        PDBG("Tolerance = "<<tolerance)
        
        getline( inNuissance, line);
        iss.clear();
        iss.str(line);
        while(iss >> doubtemp){
            delta=doubtemp;
        }
        PDBG("Delta = "<<delta)
        
        getline( inNuissance, line);
        iss.clear();
        iss.str(line);
        while(iss >> doubtemp){
            borderRange.push_back(doubtemp);
        }
        PDBG("Border Range [ "<<borderRange[0]<<", "<<borderRange[1]<<" ]")
        
        getline( inNuissance, line);
        iss.clear();
        iss.str(line);
        while(iss >> doubtemp){
            migRange.push_back(doubtemp);
        }
        PDBG("Migration Range [ "<<migRange[0]<<", "<<migRange[1]<<" ]")
        
//        getline( inNuissance, line);
//        iss.clear();
//        iss.str(line);
//        while(iss >> doubtemp){
//            aRange.push_back(doubtemp);
//        }
//        PDBG("Assortative Mating Range [ "<<aRange[0]<<", "<<aRange[1]<<" ]")

        getline( inNuissance, line);
        iss.clear();
        iss.str(line);
        while(iss >> doubtemp){
            wRange.push_back(doubtemp);
        }
        PDBG("Fitness Range [ "<<wRange[0]<<", "<<wRange[1]<<" ]")
        
        getline( inNuissance, line);
        iss.clear();
        iss.str(line);
        while(iss >> doubtemp){
            startPop=doubtemp;
        };
        PDBG("First population starts at: "<<startPop<<"km")
        
        getline( inNuissance, line);
        iss.clear();
        iss.str(line);
        while(iss >> doubtemp){
            deltaDeme=doubtemp;
        };
        PDBG("deltaDeme: "<<deltaDeme)
        
        getline(inNuissance, line);
        iss.clear();
        iss.str(line);
        while(iss >> temp){
            migDis=temp;
        };
        PDBG("migDistribution:"<<migDis);
        
        PDBG("sd for each parameter: ")
        getline( inNuissance, line);
        iss.clear();
        iss.str(line);
        while(iss >> doubtemp){
            sdparams.push_back(doubtemp);
            PDBG(doubtemp<<" ")
        };
        
    }
    
    else{
        PDBG("No 'inPhlox_nuissance.params' file, Using default pars.")
        nWorlds=1;
        wName[0]=1;
        totDemes=120;
        nZones=2;
        maxGens=50000;
        maxChange=0.00000001;
        temperature= 1.2;
        tolerance = 0.001;
        delta=1;
        temp_step=0.85;
        borderRange.push_back(0.1);borderRange.push_back(0.9);
        migRange.push_back(0);migRange.push_back(1);
        //aRange.push_back(0);aRange.push_back(1);
        wRange.push_back(0.0001);wRange.push_back(10000);
        startPop= 1.5;
        deltaDeme= 0.15;
        migDis= 1;
    };
    w.resize(nZones);
};

void Params::addUntransformedPars(vector<double> parSet){
    //cout<<"From addUntransformedPars()\n";
    //cout<<"Before : ";for(int i=0; i<9; ++i)cout<<parSet[i]<<" ";cout<<'\n';
    
    //If using logit transformation this will convert back to inverse logit and replace
    
//    border.clear(); border.resize(nWorlds);
//    for(int i=0; i<nWorlds; ++i){
//        
//        border[i]=inverseLogit(parSet[i], borderRange);
//        //PDBG("border for world "<<i<<" "<<border[i]<<" ");
//    };
//    
//	m=inverseLogit(parSet[nWorlds], migRange);
//	a=inverseLogit(parSet[nWorlds+1], aRange);
//   
//    w.clear(); w.resize(2);
//    for(int i=(nWorlds+2); i<(nWorlds+5); ++i) w[0].push_back(inverseLogit(parSet[i], wRange));
//    for(int i=(nWorlds+5); i<(nWorlds+8); ++i) w[1].push_back(inverseLogit(parSet[i], wRange));
//    
//    for(int i=0;i<nWorlds;++i){if(border[i]!=border[i]){ cerr<<"Border out of range\n";exit(1);}}
//    if(m!=m){ cerr<<"Migration out of range\n";exit(1);}
//    if(a!=a){ cerr<<"Assortative mating out of range\n";exit(1);}
//    for(int i=0;i<2;++i){for(int j=0;j<3;++j){if(w[i][j]!=w[i][j]){ cerr<<"Selection coefs out range\n";exit(1);}}}
 
    //This is for untransformed parameters (although the assortative mating is still log transformed)
    border.clear(); border.resize(nWorlds);
    for(int i=0; i<nWorlds; ++i){
        
        border[i]=parSet[i];
    };
    
	m=parSet[nWorlds+1];
    
	a=parSet[nWorlds];

    w.clear(); w.resize(2);
    for(int i=(nWorlds+2); i<(nWorlds+5); ++i) w[0].push_back(parSet[i]);
    for(int i=(nWorlds+5); i<(nWorlds+8); ++i) w[1].push_back(parSet[i]);
    
    for(int i=0;i<nWorlds;++i){if(border[i]!=border[i]){ cerr<<"Border out of range\n";exit(1);}}
    if(m!=m){ cerr<<"Migration out of range\n";exit(1);}
    //if(a!=a){ cerr<<"Assortative mating out of range\n";exit(1);}
    for(int i=0;i<2;++i){for(int j=0;j<3;++j){if(w[i][j]!=w[i][j]){ cerr<<"Selection coefs out range\n";exit(1);}}}
    


  };
vector<double> Params::VecInverseLogit(vector<double> p){
    vector<double> out; out.reserve(p.size());
    for(int i=0;i<nWorlds;++i){
        out[i]=inverseLogit(p[i],borderRange);
    };
    out[nWorlds]=inverseLogit(p[nWorlds],migRange);
    //out[nWorlds+1]=inverseLogit(p[nWorlds+1],aRange);
    for(int i=(nWorlds+1);i<(nWorlds+7);++i) out[i]=inverseLogit(p[i],wRange);
    
    return out;
};


