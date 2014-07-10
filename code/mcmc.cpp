//
//  mcmc.cpp
//  phloxMCMC
//
//  Created by Robin Hopkins on 8/8/13.
//  Copyright (c) 2013 Robin Hopkins. All rights reserved.
//
//#define MCMCDBGR
//#define SEEDDBGR
#include "mcmc.h"
#include <math.h>       /* log */

Mchain::Mchain(Params* GlobalPars) : dump("mcmc_dump.log")
{
    //this reads in the observed data used
    callcounter=0;
    global = GlobalPars;
    //dump.open("mcmc_dump.log");
    dump<<"a b1 b2 b3 b4 b5 as mig s1_db s1_dr s1_lr s2_db s2_dr s2_lr likelihood"<<std::endl;
    
    dataFromFile.resize(global->nWorlds);
    for (int i=0; i<global->nWorlds;i++){
        vector<ObsData> tempData;
        tempData= readObsFiles(global->wName[i]);
        for (int j=0; j<tempData.size();++j){
            dataFromFile[i].push_back(tempData[j]);
        };
    };
    
    seedvec=inseeds();
    totRuns=seedvec.size();
    VecDoub paramCase=startingSet();
    
    current.set=paramCase;
    current.likeli=oneRun(current.set);
    //current.likeli=oneRunPrior(current.set);
}

vector<vector<double> > Mchain::inseeds(){
    vector<vector<double> > outseeds;
    int sd_cnt=0;
    
    stringstream seeds;
    seeds<<"Rseeds.txt";
    ifstream seedfile(seeds.str().c_str());
    if(!seedfile.is_open()) cout<<"loading seeds failed failed\n";
    
    SEEDDBG("\nLoading seeds\n");
    string Loader1;
    //getline(seedfile,Loader1);
    SEEDDBG(" the seeds are: ");
    while(std::getline(seedfile,Loader1)){
                
        istringstream dataline(Loader1);

        vector<double> tempvec;
        double tempdb;

        while(dataline >> tempdb){
        tempvec.push_back(tempdb);
        SEEDDBG(tempdb<<" ");
        }
        outseeds.push_back(tempvec);
        SEEDDBG("\n");
        ++sd_cnt;
    }
    SEEDDBG("total seeds: "<<outseeds.size()<<" or "<<sd_cnt<<'\n');
    return outseeds;
}



VecDoub Mchain::startingSet(){
    VecDoub vec(global->nWorlds+8);
    int nW=global->nWorlds;
    stringstream general;
    general <<"seedSet.params";
    ifstream fin( general.str().c_str());
    if(!fin.is_open()){cerr<<"'seedSet.params' file is missing!\n"; exit(1);}
    
    string line;
    double doubtemp;
    vector<double> set;
    
    getline( fin, line);
    istringstream iss(line);
    while(iss >> doubtemp){
        set.push_back(doubtemp);
    }
    
    //This is to Logit Transform all the parameters
//    
//    for (int i=0; i<nW;++i){
//        vec[i]=Logit(set[i], global->borderRange);
//        //cerr<<"world "<<i<<" border "<<vec[i]<<"\n";
//    }
//    vec[nW]=Logit(set[nW], global->migRange);
//	vec[(nW+1)]=Logit(set[(nW+1)], global->aRange);
//    for(int i=(nW+2); i<(nW+5); ++i) vec[i]=Logit(set[i], global->wRange);
//    for(int i=(nW+5); i<(nW+8); ++i) vec[i]=Logit(set[i], global->wRange);
//    
//    cout<<"Seed parameter set\n";
//    for(int i=0; i<NPARs; ++i)   cout<<vec[i]<<" ";
//    cout<<'\n';
  
    //this is for untransformed paramters
    
    cout<<"Seed parameter set\n";
    for (int i=0;i<set.size();++i){
        vec[i]=set[i];
        cout<<vec[i]<<" ";
    }
     
    cout<<'\n';
    return vec;
}

double Mchain::oneRun(VecDoub & pars_of_interest){
    callcounter++;
    
    //if(callcounter%10==0) cout<<"|" << std::flush;
    //if(callcounter%1000==0) cout<<callcounter<< std::flush;
    //cout<<"From OneGo():\n";
    //for(int i=0; i<9; ++i) cout<<pars_of_interest[i]<<" ";cout<<endl;
    
    vector<double>likeList;
    likeList.resize(global->nWorlds);
    double likelihood=0.0;
    
    //Add the parameters of interest to the Params struct
    global->addUntransformedPars(NR_to_vec(pars_of_interest));
    
    // Create multiple instances of the World class with the parameters from which we'll run the simulation 
    for (int i=0; i<global->nWorlds; ++i){
        
        World phlox(global, i);
        
        // Run simulation!
        phlox.simulate(global, i);
        //phlox.outputPhenoF();
        
        //Get the likelihood for this parameter set
        vector<DemeOut>simDemesOut=phlox.MakesimDemeOut();
        double templike=calcTotLike(simDemesOut, dataFromFile[i], global);
        likeList[i]=templike;
        //dump<<"tran"<<global->wName[i]<<" "<<likeList[i]<<" "<<phlox.output_Gens()<<" ";
    };
    for (int i=0; i<likeList.size();++i){
        likelihood+=likeList[i];
    };
    
    //    dump<<-likelihood<<" ";
    //    vector<double> pars_for_out=global->VecInverseLogit(NR_to_vec(pars_of_interest));
    //    for(int i=0; i<NPARs; ++i)   dump<<pars_for_out[i]<<" ";
    //    dump<<endl;
    MCMCDBG("likelihood from oneRun: "<<likelihood<<'\n');
    return likelihood;
}
//This is just like oneRun but it adds a prior to the likelihood (prior function in Utilities, parameters defined above)
double Mchain::oneRunPrior(VecDoub & pars_of_interest){
    callcounter++;

    vector<double>likeList;
    likeList.resize(global->nWorlds);
    double likelihood=0.0;
    
    //Add the parameters of interest to the Params struct
    global->addUntransformedPars(NR_to_vec(pars_of_interest));
    
    // Create multiple instances of the World class with the parameters from which we'll run the simulation
    for (int i=0; i<global->nWorlds; ++i){
        
        World phlox(global, i);
        
        // Run simulation!
        phlox.simulate(global, i);
        //phlox.outputPhenoF();
        
        //Get the likelihood for this parameter set
        vector<DemeOut>simDemesOut=phlox.MakesimDemeOut();
        double templike=calcTotLike(simDemesOut, dataFromFile[i], global);
        likeList[i]=templike;
        //dump<<"tran"<<global->wName[i]<<" "<<likeList[i]<<" "<<phlox.output_Gens()<<" ";
    };
    for (int i=0; i<likeList.size();++i){
        likelihood+=likeList[i];
    };
    MCMCDBG("priors added:"<<'\n');
    for (int i=0; i<3; ++i){
        double param=global->w[1][i];
        MCMCDBG(" w: "<<param<<" ");
        double priorprob=calcPrior(param, -3, 7);
        MCMCDBG(priorprob);
        likelihood+=priorprob;
    };
    MCMCDBG('/n');
     MCMCDBG("likelihood from oneRunPrior: "<<likelihood<<'\n');
    return likelihood;
}


//This function takes the old set of parameters and makes a new set based on changing a random parameter and changing it drawn from a normal distribution with MCMCstep variance.
VecDoub Mchain::propose(vector<double>& sdpar){
    
    MCMCDBG("current set from propose function ");
    for(int j=0; j<global->nWorlds+8; ++j){
        MCMCDBG(current.set[j]<<" ");
    };
    
    VecDoub tempset = current.set;
    int paramChange = randint(global->nWorlds+1, global->nWorlds+7);    // WARNING WATCH FOR HARD CODING!!!!!!! pick number between 5-12
    double change =normalrand();
    if(paramChange < global->nWorlds){
        if (change < 0.0){
            change=global->totDemes;
            MCMCDBG(global->totDemes <<" ");
            tempset[paramChange] = tempset[paramChange]-(1/change);
            }
        else if (change > 0.0) {
            change=global->totDemes;
            MCMCDBG(global->totDemes <<" ");
            tempset[paramChange] = tempset[paramChange]+(1/change);
           }
    }
    else {tempset[paramChange] = tempset[paramChange] + sdpar[paramChange]*change;}
    
    MCMCDBG("rand int paramchange "<<paramChange<<" change "<<change<<" sd "<<sdpar[paramChange]<<" temp set "<<tempset[paramChange]<<'\n');
    return tempset;
};


VecDoub Mchain::proposeLN(vector<double>& sdpar){
    
    MCMCDBG("current set from propose function ");
    for(int j=0; j<global->nWorlds+7; ++j){
        MCMCDBG(current.set[j]<<" ");
    };
    
    VecDoub tempset = current.set;
    int paramChange = randint(0, global->nWorlds+6);    // WARNING WATCH FOR HARD CODING!!!!!!! pick number between 5-11
    //paramChange=11; //HARDCODE TO CHECK ONE PARAMETER AT A TIME!!!!!
    double change =normalrand();
    if(paramChange < global->nWorlds){
        if (change < 0.0){
            change=global->totDemes;
            MCMCDBG(global->totDemes <<" ");
            tempset[paramChange] = tempset[paramChange]-(1/change);
        }
        else if (change > 0.0) {
            change=global->totDemes;
            MCMCDBG(global->totDemes <<" ");
            tempset[paramChange] = tempset[paramChange]+(1/change);
        }
    }
    else {tempset[paramChange] = log(tempset[paramChange]) + sdpar[paramChange]*change;
        tempset[paramChange]=exp(tempset[paramChange]);}
    
    MCMCDBG("rand int paramchange "<<paramChange<<" change "<<change<<" sd "<<sdpar[paramChange]<<" temp set "<<tempset[paramChange]<<'\n');
    return tempset;
};
//Call this to propose change to selection and migration
VecDoub Mchain::multiPropose(int step){
    MCMCDBG("current set from propose function ");
    for(int j=5; j<12; ++j){
        MCMCDBG(current.set[j]<<" ");
    };
    
    VecDoub tempset = current.set;
    Doub temp;
    for(int k=0; k<seedvec[1].size(); ++k){
        temp=seedvec[step][k];
        tempset[k+ParamStart]=tempset[k+ParamStart] + temp;
        MCMCDBG("\nparam number "<<k<<" change "<<temp<< "to become "<<tempset[k+ParamStart]);
    }
    MCMCDBG('\n');
    return tempset;
}

VecDoub Mchain::multiBord(int step){
    MCMCDBG("current set from propose function ");
    for(int j=0; j<global->nWorlds+7; ++j){
        MCMCDBG(current.set[j]<<" ");
    };
    
    VecDoub tempset = current.set;
    Doub temp;
    
    //PUT IN BORDER CHECK HERE!!!!!
    for(int k=0; k<global->nWorlds; ++k){
        if (seedvec[step][k]!=0){
            seedvec[step][k]/=global->totDemes;
        }
    }
   
    for(int k=0; k<seedvec[1].size(); ++k){
        temp=seedvec[step][k];
        tempset[k]+=temp;
        MCMCDBG("\nparam number "<<k<<" change "<<temp<< " to become "<<tempset[k]);
    }
    MCMCDBG('\n');
    return tempset;
}

double Mchain::mcmcStep(){
    double accept=0;
    
    for (int i=0; i<Mlink; ++i){
        //first make the new proposed state
        VecDoub proposSet = propose(global->sdparams);
        State proposal;
        
        //check to confirm proposal set is within the biological bounds (positive)
        double negcheck=1;
        double alph;
        double rand=0;
        for(int k=0; k<proposSet.size();++k){
            
            if (proposSet[k]<0) {
                negcheck=-1;}
        }
        if (negcheck<0){
            alph=negcheck;
        }
        else {
            proposal.set=proposSet;
            proposal.likeli=oneRunPrior(proposal.set);
            MCMCDBG("current ");
            for(int j=0; j<global->nWorlds+9; ++j){
                MCMCDBG(current.set[j]<<" ");
            };
            MCMCDBG(current.likeli<<'\n');
            MCMCDBG("proposal ");
            for(int i=0; i<global->nWorlds+9; ++i){
                MCMCDBG(proposal.set[i]<<" ");
            };
            MCMCDBG(proposal.likeli<<'\n');
            
            //then decided if you're keeping the new state
            
            alph = min(1.0,exp(proposal.likeli-current.likeli));
            rand =randreal();
        };
        
        MCMCDBG("alph "<<alph<<" rand "<<rand<<"\n");
        if (rand<alph) {
            current=proposal;
            accept++;
        };
    }
    
    MCMCDBG(accept/Mlink<<'\n');
    return accept/Mlink;
};


double Mchain::mcmcMultiStep(int nsteps){
    double accept = 0;
    
    for (int i=0; i<nsteps; ++i){
        VecDoub proposSet = multiPropose(i);
        State proposal;
        
        double negcheck = 1;
        double alph;
        double rand = 0;
        for (int j=0; j<proposSet.size(); ++j){
            if (proposSet[j]<0) {
                negcheck=-1;}
        }
        if (negcheck<0){
            alph=negcheck;
        }
        else {
            proposal.set=proposSet;
            proposal.likeli=oneRun(proposal.set);
            MCMCDBG("current ");
            for(int j=global->nWorlds; j<global->nWorlds+7; ++j){
                MCMCDBG(current.set[j]<<" ");
            };
            MCMCDBG(current.likeli<<'\n');
            MCMCDBG("proposal ");
            for(int i=global->nWorlds; i<global->nWorlds+7; ++i){
                MCMCDBG(proposal.set[i]<<" ");
            };
            MCMCDBG(proposal.likeli<<'\n');
            
            //then decided if you're keeping the new state
            
            alph = min(1.0,exp(proposal.likeli-current.likeli));
            rand =randreal();
        };
        
        MCMCDBG("alph "<<alph<<" rand "<<rand<<"\n");
        if (rand<alph) {
            current=proposal;
            accept++;
        };
        printset();
    }
    cout<<"acceptance rate"<<accept/nsteps<<'\n';
    
    return (accept/nsteps);
}

double Mchain::mcmcbord(int nsteps){
    double accept = 0;
    
    for (int i=0; i<nsteps; ++i){
        VecDoub proposSet = multiBord(i);
        State proposal;
        
        double negcheck = 1;
        double alph;
        double rand = 0;
        for (int j=global->nWorlds; j<proposSet.size(); ++j){
            if (proposSet[j]<0) {
                negcheck=-1;}
        }
       
        if (negcheck<0){
            alph=negcheck;
        }
              
        else {
            proposal.set=proposSet;
            proposal.likeli=oneRun(proposal.set);
            MCMCDBG("current ");
            for(int j=0; j<(global->nWorlds+7); ++j){
                MCMCDBG(current.set[j]<<" ");
            };
            MCMCDBG(current.likeli<<'\n');
            MCMCDBG("proposal ");
            for(int i=0; i<(global->nWorlds+7); ++i){
                MCMCDBG(proposal.set[i]<<" ");
            };
            MCMCDBG(proposal.likeli<<'\n');
            
            //then decided if you're keeping the new state
            
            alph = min(1.0,exp(proposal.likeli-current.likeli));
            rand =randreal();
        };
        
        MCMCDBG("alph "<<alph<<" rand "<<rand<<"\n");
        if (rand<alph) {
            current=proposal;
            accept++;
        };
        printset();
    }
    cout<<"acceptance rate"<<accept/nsteps<<'\n';
    
    return (accept/nsteps);
}




void Mchain::printset(){
    //This is printset if logit transformation was done
    
//    vector<double> unTranSet= global->VecInverseLogit(NR_to_vec(current.set));
//    for(int i=5; i<13; ++i){
//        dump<<unTranSet[i]<<" ";
//    };
    vector<double> outparams=NR_to_vec(current.set);
    
    for(int i=0; i<(global->nWorlds+8); ++i){
        
        dump<<outparams[i]<<" ";
    };

    dump<<current.likeli<<std::endl;
    
    
};


