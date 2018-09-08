//
//  population.hpp
//  sandglass
//
//  Created by Nate Layman on 1/4/17.
//  Copyright © 2017 Nate Layman. All rights reserved.
//

#ifndef population_h
#define population_h

#include "individual.hpp"

#include <random>
#include <valarray>
#include <unordered_set>

class population {
    
public:
    
    std::string popName;
    unsigned int ploidy,nloci, sCounter, domLevels, model, SAlleleNumber;
    double inbreeding_depression,SModifierFreq,qbar,half_Max,fitnessSum,Mu,s,h,probSelf,effort;
    std::valarray<double> genotypeFreqs;
    std::mt19937 p_rng;
    std::uniform_real_distribution<> realDist;
    std::uniform_int_distribution<> FDR_map;
    std::vector<individual> parents, offspring;
    std::unordered_set<int> SAllele_set;
    
    population(std::mt19937& rng,
               unsigned int _ploidy,
               unsigned int _N,
               unsigned int _nloci,
               unsigned int _sCounter,
               unsigned int _domLevels,
               unsigned int _model,
               double _s,
               double _h,
               double _Mu,
               double _probSelf); // Constructor
    
    population( const population &obj); // Copy constructor
      
    void growPop(double r, double K, bool diffusion);
    individual getGamete(individual gamete, bool unreduced);
    bool cross(bool wgd, unsigned int progenyID);
    unsigned int chooseParent(void);
    bool reproduce(bool wgd); // Have this return inbreeding depression after reproduction
    void adjustPopSize(int popSize);
    int getGamete(bool wgd, bool FDR, unsigned int mutations);
    void mutate(bool burn, bool reversal);
    void printParents(void);
    void printOffspring(void);
    bool competePhenotypes( std::vector<std::vector<int>>& mother, std::vector<std::vector<int>>& father );
    
};

#endif // POPULATION_H

