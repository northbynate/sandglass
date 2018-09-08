//
//  locus.cpp
//  sandglass
//
//  Created by Nate Layman on 1/4/17.
//  Copyright © 2017 Nate Layman. All rights reserved.
//

#include "locus.hpp"
#include <iostream>

// Constructor
locus::locus(double* _s, double* _h) :
s(_s),
h(_h),
fitness(1),
mutations(0) {}
 
void locus::mutate(double half_Max, int ploidy)
{
    mutations++;
    getFitness(half_Max, ploidy);
}

double locus::getFitness( double half_Max, int ploidy )
{
    // Dominance:
    double X = (double)mutations / ploidy;
    
    double hsubx = 0;
    if(*h && X && X < 1) hsubx = 1 / ( 1 + half_Max * ( 1 - X ) / X );
    if(X >= 1) hsubx = 1;
     
    // Fitness:
    fitness = 1 - hsubx * *s;
    
    if(ploidy > 2)
    {
    //std::cout << "ploidy: " << ploidy << " mutations: " << mutations << " X: " << X << " h: " << hsubx << " s: " << *s << " fitness: " << fitness << std::endl;
    }
        
    return fitness;
}

