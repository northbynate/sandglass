//
//  locus.hpp
//  sandglass
//
//  Created by Nate Layman on 1/4/17.
//  Copyright © 2017 Nate Layman. All rights reserved.
//

#ifndef locus_hpp
#define locus_hpp

class locus {
      
public:
    void mutate(double half_Max, int ploidy);
    
    double fitness, *s, *h;
    unsigned int mutations;
    
    locus(double* _s, double* _h); //Constructor

    double getFitness(double half_Max, int ploidy);
    
};

#endif

