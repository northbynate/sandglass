//
//  main.cpp
//  sandglass
//
//  Created by Nate Layman and Jeremiah Busch on 1/4/17.
//  Copyright © 2017 Nate Layman and Jeremiah Busch. All rights reserved.
//
// 

/* model options:
 0: random mating
 1: SSI - codom
 2: SSI - domcod: domcod model.  Pollen has dominance, pistil has codominance.
 3: SSI - dom
 4: GSI -
 All of these SI models conform to fecundity selection with t=1 trial (Vekemans et al. 1998)
 The models below are not observed naturally, but were coded for completeness:
 5: GSI - domcod: domcod model.  Pollen has dominance, pistil has codominance. =codom in diploids.
 6: GSI - dom  */

#include <omp.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include <sstream>
#include <fstream>
#include <cstring>
#include <sys/stat.h>
#include <sys/file.h>
#include <unistd.h>
#include "population.hpp"

double s = 0.1;
double h = 0.2; // Dominance. 0 = complete recessivity.
double Mu = 1;  // Mutation rate per DIPLOID genome.
double probSelf = 0;
double lambda = 1.05;
unsigned int basePloidy = 2;
unsigned int nloci = 1024;
unsigned int mainland_N = 1000;
unsigned int postBottleneck_N = 5;
unsigned int postGenerations = 10; // Max generations on island
unsigned int iters = 1000; // independent runs using the same parameters
unsigned int window = 35; // Anything over 50 takes forever
unsigned int domLevels = 1;
unsigned int sRichness = 40; // Start at 30 so we equilibrate faster (Vekemans et al. 1998)
unsigned int model = 0;

bool brief = false;   // brief trumps logistic growth.
bool logistic_growth = true;
bool diffusion=false;    // Adds a random element to logistic growth. Demographic variation

bool hard_selection=false;
bool wgd=false;
bool two_way=false; // Permits two way mutation at the M loci

std::string filename = "results.csv";

void writeOutput(std::stringstream &results, population &popN, int iterations, int generations);
void collapsePopN(population &popN);
void displaySettings();
bool parseOptions(int argc, char* argv[]);
inline bool fileExists (const std::string& name);

int main(int argc, char* argv[])
{
    std::random_device rd{};
    
    //Unsigned in case of overflow
    uint32_t seed = rd();
    
    if(parseOptions(argc, argv)) return 5;
    
    //Sanitize some of the inputs
    nloci += 1; // Account for the M-locus. It is always the first one in the vector (locus 0).
    if(!domLevels) domLevels = 1;
    if(!sRichness) sRichness = 1;
    if (basePloidy % 2) basePloidy += 1; // Make sure basePloidy is an even number
    
    std::ofstream resultsFile;
    
    resultsFile.open(filename , std::ios::trunc);
    resultsFile << "pop,model,lambda,iteration,generation,ploidy,s,h,Mu,N,N.bottle,nloci,probSelf,ID,qbar,SmF,Srich,effort,wbar";
    for(int cntr = 0; cntr < (basePloidy * (1 + wgd) + 1); ++cntr ) { resultsFile << "," << cntr; }
    resultsFile << std::endl;

    resultsFile.precision(10);
    
    double start_time = omp_get_wtime();

    #pragma omp parallel
    {
		std::stringstream buff;
		
		#pragma omp for schedule(dynamic)
		for (int iterations = 0; iterations < iters; ++iterations)
		{
            
            // A note about random number generation:
            // seed is a 32 bit integer generated using
            // a random_device. Mersenne Twister is guarenteed
            // to produce a different series of random numbers
            // when fed unique seeds - even if those seeds are
            // in series. By adding the iteration number to seed
            // we avoid the birthday problem and insure that all
            // iterations get unique seeds as well as ensuring
            // that the code won't produce the same results when run
            // twice. We initialize gen inside this loop so that each
            // iteration will have it's own random number generator
            // to ensure that different threads don't have to wait
            // to access a single instance of gen.
            std::mt19937 gen(seed + iterations);

            double growthRate = log(lambda);
		
			// Form mainland population
			population mainland(gen, basePloidy,mainland_N,nloci,sRichness,domLevels,model,s,h,Mu,probSelf);
		
			bool equilibrium = false;
			bool mutateS = false;
			double previous = 0;
			double current;
			unsigned int oscillations = 0;
			unsigned int generations;
		
			for(generations = 0; !equilibrium; ++generations)
			{
				mainland.mutate(mutateS, two_way ); // Don't mutate the selfing modifier until burn in has passed.
			
				// If true, returned population went extinct due to lack of compatible mates.
				if( mainland.reproduce(false) ) { collapsePopN(mainland); break; }
			
				current = mainland.inbreeding_depression; // Wait for ID equilibrium before introducing the modifier
			
				// even oscillations happen when the sign has moved from negative (or zero) to positive.
				if( !(oscillations % 2) && (current - previous) > 0) oscillations ++;
				// odd oscillations happen when the sign has moved from positive to negative.
				else if( (oscillations % 2) && (current - previous) < 0) oscillations ++;
			
				writeOutput(buff, mainland, iterations, generations);
						
				if( !(generations % window) ) { // how wide is the oscillation window?
				
					// To figure out if a population is at equilibrium check
					// the number of oscillations inside a window of time.
					// An oscillation is when the sign of the difference between
					// the first generation in the window and any subsequent generation
					// changes.
					if(oscillations >= window/2)
					{
						if(!mutateS) mutateS = true;
						else equilibrium = true;
		                        } else {
						oscillations = 0;
						previous = current;
					}
				}
			
			} // end generations for loop
		
			// Mutate the mainland one last time.
			mainland.mutate(true, two_way);
		
			// Copy the mainland population into a new population which will experience a bottleneck.
			population island(mainland);
			island.popName = "isl";
		
			// Form the bottleneck. This puts a limit on the number of offspring that can be made in the next reproduce().
			island.adjustPopSize(postBottleneck_N);
			island.reproduce(false);
		
			population wgd_island(mainland);
			if(wgd) {
				wgd_island.popName = "wgd_isl";
				wgd_island.adjustPopSize(postBottleneck_N);
				wgd_island.reproduce(wgd);
			}
		
			writeOutput(buff, mainland, iterations, generations);
		
			// Run post bottleneck populations
			for(generations = 0; generations < postGenerations; ++generations)
			{
				if(island.parents.size())
				{
					if(island.offspring.size())
					{
				 		island.mutate(true, two_way); // First true allows mutations to selfing modifier.
				 		if(island.reproduce(false)) collapsePopN(island); // False because this popN doesn't experience WGD.
				 		if( brief && island.parents.size() < mainland_N) island.adjustPopSize(mainland_N);
				 		else if (logistic_growth) island.growPop(growthRate,mainland_N,diffusion);
				 	} else collapsePopN(island);
				 	writeOutput(buff, island, iterations, generations);
				}
			
				if(wgd && wgd_island.parents.size())
				{
					if(wgd_island.offspring.size())
					{
						wgd_island.mutate(true, two_way); // First true allows mutations to selfing modifier.
						if(wgd_island.reproduce(false)) collapsePopN(wgd_island); // False because no bottleneck is going on. Genetic extinction.
						if( brief && wgd_island.parents.size() < mainland_N) wgd_island.adjustPopSize(mainland_N);
						else if (logistic_growth) wgd_island.growPop(growthRate,mainland_N,diffusion);
					} else collapsePopN(wgd_island);
					writeOutput(buff, wgd_island, iterations, generations);
				}
			
			} // Close generations for loop

		#pragma omp critical
                {
                    resultsFile << buff.rdbuf();
                }

	    } // Close iterations loop
       
    } // Close omp parallel  
	
    resultsFile.close();
	
    std::ifstream in(filename, std::ifstream::ate | std::ifstream::binary);
    auto size = in.tellg();
	
    double time = (omp_get_wtime() - start_time) / 60; // in minutes
    double rate = size / (time * 1000000); // 1,000,000 converts from bytes to megabytes
    
    std::cout << "filename: " << filename << " time: " << time << " size: " << size << " rate: " << rate << std::endl;
    return 0;
}

void writeOutput(std::stringstream &results, population &popN, int iterations, int generations)
{
    
    results << popN.popName << ","
    << model << ","
    << lambda << ","
    << iterations << ","
    << generations << ","
    << popN.ploidy << ","
    << popN.s << ","
    << popN.h << ","
    << popN.Mu << ","
    << popN.parents.size() << ","
    << postBottleneck_N << ","
    << popN.nloci-1 << ","
    << popN.probSelf << ","
    << popN.inbreeding_depression << ","
    << popN.qbar << ","
    << popN.SModifierFreq << ","
    << popN.SAlleleNumber << ","
    << popN.effort;
    
    double wbar = popN.parents.size() ? (popN.fitnessSum / popN.parents.size()) : 0;
    results << "," << wbar;
    
    for(unsigned int geno = 0; geno < popN.genotypeFreqs.size(); ++geno) results << "," << popN.genotypeFreqs[geno];
    results << "\n";
    
    // Uncomment below for debuging. It will circumvent writing to the file in lieu of displaying results via std::cout.
    // std::cout << results.rdbuf();
}

void collapsePopN(population &popN)
{
    popN.parents.clear();
    popN.offspring.clear();
    popN.inbreeding_depression = 0;
    popN.qbar = 0;
    popN.SModifierFreq = 0;
    popN.SAlleleNumber = 0;
}

void displaySettings()
{
    std::cout <<  std::endl <<
    "            Sandglass by Nate Layman and Jeremiah Busch" <<  std::endl << std::endl <<
    " Settings:" << std::endl << std::endl <<
    "   -?            Displays this help file" << std::endl <<
    "   -brief        Enables instant population size recovery after the bottleneck. Default: false" << std::endl <<
    "   -diffusion    Introduces a random element to logistic growth. Default: true" << std::endl <<
    "   -dom [#]      Number of dominance levels. Default: 1" << std::endl <<
    "   -f filename   Output file name. Default: results.csv " << std::endl <<
    "   -g [#]        Post-bottleneck generations. Default: 250" << std::endl <<
    "   -lambda [#.#]     Net reproductive rate during popN recovery. Default: 1.02" << std::endl <<
    "   -h [0.#]      Dominance coefficient. Default: 0.2" << std::endl <<
    "   -i [#]        Iterations. Default: 100" << std::endl <<
    "   -l [#]        Number of ID loci. Default: 1024" << std::endl <<
    "   -logic        Enables logistic growth after the bottleneck. Default: true" << std::endl <<
    "   -model [#]    0 Random mating. This disables the S-locus" << std::endl <<
    "                 1 SSI codom" << std::endl <<
    "                 2 SSI domcod" << std::endl <<
    "                 3 SSI dom" << std::endl <<
    "                 4 GSI codom" << std::endl <<
    "                 5 GSI domcod" << std::endl <<
    "                 6 GSI dom" << std::endl <<
    "                 default: 1" << std::endl <<
    "   -Mu [#.#]     Diploid genome mutation rate. Default: 1" << std::endl <<
    "   -N [#]        Mainland starting population size. Default: 1000" << std::endl <<
    "   -neo          Enable polyploidization during the bottleneck. Default: false" << std::endl <<
    "   -p [#]        Mainland ploidy. Default: 2 " << std::endl <<
    "   -post [#]     Post bottleneck population size. Default: 20"
    "   -s [0.#]      Selection coefficient. Default: 0.1" << std::endl <<
    "   -S [0.#]      Modifier selfing rate. 0 disables the selfing modifier. Default: 0" << std::endl <<
    "   -Srich [#]    Starting number of S-alleles. Default: 5" << std::endl <<
    "   -w [#]        Stability window size in generations. Default: 20" << std::endl << std::endl <<
    "   Notes: To enable the selfing modifier set -S to a nonzero value and select model 0." << std::endl <<
    "          If brief and logic are false population will remain small following bottleneck." << std::endl << std::endl;
}

bool parseOptions(int argc, char* argv[])
{
    for(int a = 1; a < argc; a++)
    {
        if(strcmp(argv[a],"-p") == 0 && strncmp(argv[a+1],"-",1)) { basePloidy = atoi(argv[a+1]); a++; }
        else if(strcmp(argv[a],"-f") == 0 && (a+1) <= argc && strncmp(argv[a+1],"-",1)) { filename = argv[a+1]; a++; }
        else if(strcmp(argv[a],"-h") == 0  && strncmp(argv[a+1],"-",1)) { h = atof(argv[a+1]); a++; }
        else if(strcmp(argv[a],"-S") == 0  && strncmp(argv[a+1],"-",1)) { probSelf = atof(argv[a+1]); a++;  }
        else if(strcmp(argv[a],"-s") == 0  && strncmp(argv[a+1],"-",1)) { s = atof(argv[a+1]); a++; }
        else if(strcmp(argv[a],"-Mu") == 0  && strncmp(argv[a+1],"-",1)) { Mu = atof(argv[a+1]); a++; }
        else if(strcmp(argv[a],"-N") == 0 && strncmp(argv[a+1],"-",1)) { mainland_N = atoi(argv[a+1]); a++; }
        else if(strcmp(argv[a],"-post") == 0 && strncmp(argv[a+1],"-",1)) { postBottleneck_N = atoi(argv[a+1]); a++; }
        else if(strcmp(argv[a],"-w") == 0 && strncmp(argv[a+1],"-",1)) { window = atoi(argv[a+1]); a++; }
        else if(strcmp(argv[a],"-i") == 0 && strncmp(argv[a+1],"-",1)) { iters = atoi(argv[a+1]); a++; }
        else if(strcmp(argv[a],"-g") == 0 && strncmp(argv[a+1],"-",1)) { postGenerations = atoi(argv[a+1]); a++; }
        else if(strcmp(argv[a],"-l") == 0 && strncmp(argv[a+1],"-",1)) { nloci = atoi(argv[a+1]); a++; }
        else if(strcmp(argv[a],"-brief") == 0) brief = true;
	else if(strcmp(argv[a],"-nbrief") == 0) brief = false;
        else if(strcmp(argv[a],"-?") == 0) {displaySettings(); return true;}
        else if(strcmp(argv[a],"-neo") == 0) wgd = true;
        else if(strcmp(argv[a],"-logic") == 0) logistic_growth = true;
	else if(strcmp(argv[a],"-nlogic") ==0) logistic_growth = false;
        else if(strcmp(argv[a],"-diffusion") == 0) diffusion = true;
        else if(strcmp(argv[a],"-ndiffusion") ==0) diffusion = false;
        else if(strcmp(argv[a],"-model") == 0 && strncmp(argv[a+1],"-",1)) { model = atoi(argv[a+1]); a++; }
        else if(strcmp(argv[a],"-lambda") == 0  && strncmp(argv[a+1],"-",1)) { lambda = atof(argv[a+1]); a++; }
        else if(strcmp(argv[a],"-dom") == 0 && strncmp(argv[a+1],"-",1)) { domLevels = atoi(argv[a+1]); a++; }
        else if(strcmp(argv[a],"-rich") == 0 && strncmp(argv[a+1],"-",1)) { sRichness = atoi(argv[a+1]); a++; }
        else if(strcmp(argv[a],"-2way") == 0) two_way = true;
        else if(strcmp(argv[a],"-1way") == 0) two_way = false;
    }
    return false;
}

inline bool fileExists (const std::string& name) {
    std::ifstream f(name.c_str());
    return f.good();
}


/*
 
 Fixed / working:
 mating rules (random, co-dom, domcod, dom)
 SSI mutation
 SSI dominance mutation
 Implement a way to calcuate S-allele richness each generation. Use unordered set.
 Double check how S-locus inheritence works
 Finish implementing logistic growth. Add growth rate variable.
 Add polyploidization bool to save time.
 Figure out what happens when populations collapse in postbottleneck populations. Sometimes pop size doesn't grow? Offspring size set to zero.
 Implement file output
 Help file
 Double check srichness is going down following bottleneck.
 SmF generation bullshit
 Check against Vekeman's
 Fixed unordered set error
 Continue on to next iteration if mainland pop goes extinct
 Don't continue reproducing a population if it goes extinct in post-bottleneck
 Don't continue reproducting a population if it reaches carrying capacity.
 Double check - do I choose both a new mother and a new father every time? Is this appropriate?
 slocus:
 -mutations
 -initalize
 -error check
 Save some estimation of the ammount of pollen wasted in incompatible crosses.
 nan number error during island population collapse.
 Removed condition requiring s-modifier to survive bottleneck. Bottleneck events used to be detected by parents > offspring. With diffusion that doesn't work.
 Make sure gen following extinction has N=0 then stops.
 Implement mpirun protections against file writing race conditions
 Fix effort
 Something is wrong with the genotype frequency counter.
 qbar works - verified against genotypeFreqs
 Fitness doesn't seem to vary as much as it should.
 Seems to require a pretty high growth rate just to offset selection.
 main,1,0.0953102,1.1,3,34,2,0.1,0,1,0,5,1024,0,-inf,nan,nan,0,inf,nan,nan,nan
 effort, qbar, ID all need to be stored in offspring as local variable. Needed to deal with dynamic offspring number
 Make sure s-locus is updated during wgd
 SOMETHING IS WRONG WITH BOTTLENECK. TOO MUCH ID IS LOST. CHECK CLEARING PARENTS ECT...
 
 Left to do:
 mean fitness sometimes goes over 1. Figure that out.
 */
