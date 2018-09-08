//
//  individual.cpp
//  sandglass
//
//  Created by Nate Layman on 1/4/17.
//  Copyright © 2017 Nate Layman. All rights reserved.
//

#include "locus.hpp"
#include "individual.hpp"
#include <iostream>
#include <random>
 
individual::individual( unsigned int* _ploidy,
                       unsigned int* _nloci, double* _s, double* _h ) :
ploidy(_ploidy),
nloci(_nloci),
mutations(0),
fitness(1),
s(_s),
h(_h),
ID(0),
genome( *_nloci,{s, h} ),
slocus( *ploidy, std::vector<int>(2)) {} // The mutation variable at this level tracks the total number of deleterious mutations in the indiviual.

void individual::mutate(unsigned int targetLocus, double half_Max, bool reversal, std::mt19937& gen)
{
    // Note: This mutation function assumes an equal probability of mutation for both selfing modifiers
    // and inbreeding alleles. Both get mutated using the same function.

    std::uniform_real_distribution<> realDist(0,1);
    // The mutation variable at this level refers to mutations responsible for inbreeding depression. It doesn't include the S-locus.
    if(targetLocus)
    {
        mutations++;
        
        // Go ahead with mutating the target locus.
        genome[targetLocus].mutate(half_Max, *ploidy);
        
    } else // If we're talking about the M-locus then back mutations may be allowed.
    {
        double ratio = (double)genome[targetLocus].mutations / * ploidy;
        double bar = realDist(gen);
        
        if( bar > ratio) genome[targetLocus].mutations++;
        else if( genome[targetLocus].mutations && reversal ) genome[targetLocus].mutations--;
    }
}

int individual::getGamete(bool unreduced, unsigned int locus, std::mt19937& gen)
{
    unsigned int result = 0;
    unsigned int genotype = genome[locus].mutations;

    std::bernoulli_distribution fdr(0.5);
    
    bool FDR_SDR = fdr(gen);
    
    if(genotype)
    {
        // First figure out if gamete needs to be unreduced. If so then flip FDR SDR coin
        int trials = *ploidy / (2 - unreduced * FDR_SDR);
        if(genotype >= *ploidy) result = trials;
        else if(genotype == 0) result = 0;
        else
        {
            // Unordered sample without replacment.
            for( int c = 0; c < trials; ++c ) // Bernoulli trials with shifting p
            {
                double PofA = (double)(genotype-result)/(*ploidy-c);
                std::bernoulli_distribution bernie(PofA);
                result += bernie(gen);
            }
        }
        if (unreduced && !FDR_SDR) result *= 2; // SDR : FDR
    }
    
   // if( genotype==1 ) std::cout << result << std::endl;
    
    return result;
}
