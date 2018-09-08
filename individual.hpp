//
//  individual.hpp
//  sandglass
//
//  Created by Nate Layman on 1/4/17.
//  Copyright © 2017 Nate Layman. All rights reserved.
//
 
#ifndef individual_h
#define individual_h

#include "locus.hpp"
#include <vector>
#include <random>
#include <valarray>
 
class individual {
    
public:
    
    unsigned int* ploidy, * nloci;
    unsigned long mutations;
    double fitness, *s, *h, ID;
    std::vector<locus> genome;
    std::vector<std::vector<int>> slocus;
    std::valarray<double> genotypeCounts;
    
    void mutate( unsigned int targetLocus, double half_Max, bool reversal,std::mt19937& gen );
    void mutateS_ID( unsigned int randomChromosome, unsigned int sRichness );
    void mutateS_Dom( unsigned int randomChromosome, unsigned int dom );
    int getGamete(bool unreduced, unsigned int locus, std::mt19937& gen);
    void printGenome(void);
    
    individual( unsigned int* _ploidy, unsigned int* _nloci, double* _s, double* _h ); //Constructor
};

#endif
