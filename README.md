PhloxSelection
==============

Population genetic model estimating reinforcing selection 

README for C++ code used in Hopkins et al. “Strong Reinforcing Selection in a Texas Wildflower”

About this code:
When compiled, this program will perform a Metropolis-Hastings MCMC routine to explore the likelihood surface of the data given a model with certain parameters. This program requires the Boost C++ library in order to properly compile. Please refer to the manuscript for more details on what this program does.

We have included example input files and parameter files. 

The initial genotype frequencies for starting the simulation are calculated from the “init_Phlox_IF#.txt” files.  These files contain phenotype frequencies for each of the observed populations. There is one file for each transect being simulated (in the manuscript there are 5).

The input files “inPhlox_Phen_T#.data” contain the observed phenotype counts across the transects and the distance each population is from the start of the transect. The input files “inPhlox_Fam_T#.data” contain the observed offspring counts from each maternal plant.  These files contain the population, the maternal phenotype, and the counts of offspring flower colors. 

 
The nuisance parameters in the model are the file called “inPhlox_nuissance.params”. This file contains:

The number of transects modeled
The transect names (as numbers)
The number of zones of selection
An option to stop the run after a maximum number of generations (not implemented in default settings)
An option to change the maximum change in allele frequency allowed before equilibrium is declared
An option to set the ranges of the parameters being estimated (not implemented in default settings)
The distance along a transect at which the first sampled population is placed
The distance between demes on the modeled transect
Migration kernel used
Variance for each parameter from which the MCMC algorithm samples. 

The order of the variances is:
Borders of transect 1-5, assortative mating parameter (set to zero when assume no assortative mating), migration parameter, relative fitness in the first selection zone for phenotype 2,3,4 (corresponding to dark-blue, dark-red, light-red in the manuscript), relative fitness in the second selection zone for phenotype 2,3,4 (corresponding to dark-blue, dark-red, light-red in the manuscript). 

The “seedSet.params” file contains the starting parameters. They are in the same order as the variances and correspond to the following values in the manuscript:
Border transect 1, Border transect 2, Border transect 3, Border transect 4, Border transect 5, Assortative mating, Migration, Relative fitness of Dark-blue in allopatry, relative fitness of dark-red in allopatry, relative fitness of light-red in allopatry, relative fitness of dark-blue in sympatry, relative fitness of dark-red in sympatry, relative fitness of light-red in sympatry.  
 
The MCMC algorithm starts with the parameter values indicated in the “seedSet.params” file and samples a single parameter at at time using the variances indicated in the “inPhlox_nuissance.params” file.  

The output file for this program, “mcmc_dump.log” contains the parameter values sampled by the MCMC algorithm.  
a - acceptance (1) or rejection (0) of new proposed parameters
b1 - border for transect 1
b2 - border for transect 2
b3 - border for transect 3
b4 - border for transect 4
b5 - border for transect 5
as - assortative mating 
mig - migration variance
s1_db - relative fitness in zone 2 for dark-blue
s1_dr - relative fitness in zone 2 for dark-red
s1_lr - relative fitness in zone 2 for light-red
s2_db - relative fitness in zone 2 for dark-blue
s2_dr - relative fitness in zone 2 for dark-red
s2_lr - relative fitness in zone 2 for light-red
likelihood - log likelihood of model fit to data

