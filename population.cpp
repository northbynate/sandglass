//
//  population.cpp
//  sandglass
//
//  Created by Nate Layman on 1/4/17.
//  Copyright © 2017 Nate Layman. All rights reserved.
//
 
#include "population.hpp"
#include "individual.hpp"
#include <algorithm>
#include <iostream>
#include <random>
#include <valarray>

// Constructor
population::population(std::mt19937& rng,
                       unsigned int _ploidy,
                       unsigned int _N,
                       unsigned int _nloci,
                       unsigned int _sCounter,
                       unsigned int _domLevels,
                       unsigned int _model,
                       double _s,
                       double _h,
                       double _Mu,
                       double _probSelf) :
p_rng(rng),
popName( "main" ),
ploidy( _ploidy ),
nloci( _nloci ),
sCounter(_sCounter),
domLevels(_domLevels),
model(_model),
SAlleleNumber(0),
inbreeding_depression(0),
SModifierFreq(0),
qbar(0),
half_Max(0),
fitnessSum(_N),
Mu(_Mu),
s( _s ),
h( _h ),
probSelf ( _probSelf ),
genotypeFreqs(ploidy+1),
realDist(0,1),
FDR_map(0, 63),
parents( _N,{&ploidy, &nloci, &s, &h} ),
offspring(parents)
{
    std::uniform_int_distribution<> dom(0, _domLevels - 1);
    std::uniform_int_distribution<> id(1, _sCounter);
    
    for( unsigned int i = 0; i < parents.size() ; ++i )
    {
        for( unsigned int p = 0; p < ploidy; ++p )
        {
            parents[i].slocus[p][0] = dom(p_rng);
            parents[i].slocus[p][1] = id(p_rng);
        }
    }
    
    offspring = parents;
    
    // Translate dominance coefficient into halfmax value. This only needs to be done once when the
    // population is created.
    
    half_Max = 0;
    if(h>0.000001)
    {
        if(h < 0.5)
        {
            half_Max = 1 / ( 2 * h );
            double y = log((1/h)-1) / log(half_Max);
            half_Max = pow(half_Max, y);
        }
        else if (h > 0.499999 && h < 0.500001) half_Max = 1;
        else if ( h > 0.5 && h < 1 )
        {
            half_Max = 2 - 2 * h;
            double y = log((1/h)-1) / log(half_Max);
            half_Max = pow(half_Max, y);
        }
    }
}

population::population( const population &obj):
p_rng( obj.p_rng ),
popName( obj.popName ),
ploidy( obj.ploidy ),
nloci( obj.nloci ),
sCounter(obj.sCounter),
domLevels(obj.domLevels),
model(obj.model),
SAlleleNumber( obj.SAlleleNumber ),
inbreeding_depression(obj.inbreeding_depression),
SModifierFreq(obj.SModifierFreq),
qbar(obj.qbar),
half_Max( obj.half_Max ),
fitnessSum(obj.fitnessSum),
Mu( obj.Mu ),
s( obj.s ),
h( obj.h ),
probSelf ( obj.probSelf ),
genotypeFreqs( obj.genotypeFreqs),
realDist(obj.realDist),
FDR_map( obj.FDR_map ),
parents( obj.parents ),
offspring ( parents.size(),{&ploidy, &nloci, &s, &h})
{}

// Member functions
bool population::reproduce(bool wgd)
{
    effort = 0;
    SAllele_set.clear();
    SAlleleNumber = 0;
    SModifierFreq = 0;
    inbreeding_depression = 0;
    genotypeFreqs = 0; // Note: genotypeFreqs is an ARRAY of values. See: std::valarray for use.
    qbar = 0;

    if(offspring.size())
    {
        if(wgd)
        {
            genotypeFreqs.resize(((ploidy*2)+1),0);
            for(unsigned int i = 0; i < offspring.size(); ++i) offspring[i].slocus.clear();
            ploidy *= 2;
        }
        
        // Initiate matings.
        for(unsigned int i = 0; i < offspring.size(); ++i)
        {
            offspring[i].ID = 0;
            if(cross(wgd, i)) { effort = 1000; return true; } // If popN goes extinct return true.
            SModifierFreq += offspring[i].genome[0].mutations;
            qbar += offspring[i].mutations;
            inbreeding_depression += offspring[i].ID;
        }
        
        SModifierFreq /= (offspring.size() * ploidy); //Gives average frequency of mutants at the S-locus.
        qbar /= (offspring.size() * (nloci-1) * ploidy); //( nloci - 1 ) we don't care about the selfing modifier or the S-locus here.
        inbreeding_depression /= offspring.size();
        effort /= offspring.size();
        SAlleleNumber = (int)SAllele_set.size();
        parents = std::vector<individual>(offspring);
    } else return true;
    
    // Need to figure out fitnessSum for the new generation to determine wbar
    fitnessSum = 0;
    for(unsigned int i = 0; i < parents.size(); ++i) fitnessSum += parents[i].fitness;
    return false; // false means population didn't go extinct.
}

bool population::cross(bool wgd, unsigned int progenyID)
{
    // Randomly choose a mother using relative fitness as a selection probability
    // and get two gametes from her. One ovule and one selfed pollen grain.
    // Do the whole process twice.
    
    /* fill the gamete vector with random parents. 0 and 1 are the same.
     Should only need at most 5 parents. Will need all 6 gametes though.
     0 = IDself ovule,
     1 = IDself pollen
     2 = IDout ovule,
     3 = IDout pollen,
     4 = offspring ovule,
     5 = offspring pollen
     */
    
    std::vector<unsigned int> gametes(6);
    for(int i = 0; i < 6; ++i) gametes[i] = chooseParent();
    gametes[1] = gametes[0];
    
    // 1. Identify parents. 2. Test for SC. 3. Test for cross-compatibility. If not-compatible go back to 1 until compatible.
    std::vector<std::vector<int>> pollen_S;
    std::vector<std::vector<int>> stigma_S;
    
    unsigned int tries = 0;
    do
    {
        if(tries >= 1000) return true; // Population went extinct due to lack of compatible mates.
        
        tries++;
        gametes[4] = chooseParent();
        gametes[5] = chooseParent();
        stigma_S = parents[gametes[4]].slocus;
        
        bool self = parents[gametes[4]].genome[0].mutations && (realDist(p_rng) < probSelf);
        pollen_S = self ? stigma_S : parents[gametes[5]].slocus;
        std::shuffle( std::begin(pollen_S), std::end(pollen_S), p_rng);
        
        if(model > 3) //SDR
        {
            if(wgd && (realDist(p_rng)>0.5))
            {
                pollen_S.resize(pollen_S.size()/2);
                pollen_S.insert(std::end(pollen_S), std::begin(pollen_S), std::end(pollen_S));
            } else pollen_S.resize(pollen_S.size()/(2-wgd)); //FDR
        }
        
        if(self) {
	    gametes[5] = gametes[4];
            break;
        }
        
    } while (!competePhenotypes(stigma_S, pollen_S));
    
    effort += tries;

    // Now let's get the genotye of the mother at the S-locus
    // We didn't need it before because the stigma is made of
    // maternal tissue this will always be = ploidy even during WGD.
    // We need it now to know the S-genotype of the female gametophyte.
    std::shuffle( std::begin(stigma_S), std::end(stigma_S), p_rng );
    if(wgd && (realDist(p_rng)>0.5))
    {
        stigma_S.resize(stigma_S.size()/2);
        stigma_S.insert(std::end(stigma_S), std::begin(stigma_S), std::end(stigma_S));
    } else stigma_S.resize(stigma_S.size()/(2-wgd)); //FDR

    // If GSI we've already got a gamete from the father at the s-locus. If not get it.
    if(model<=3)
    {
        if(wgd && (realDist(p_rng)>0.5))
        {
            pollen_S.resize(pollen_S.size()/2);
            pollen_S.insert(std::end(pollen_S), std::begin(pollen_S), std::end(pollen_S));
        } else pollen_S.resize(pollen_S.size()/(2-wgd)); //FDR
    }

    // If we've made it this far mating was successful. Start by figuring out offspring S-locus.
    stigma_S.insert(std::end(stigma_S),std::begin(pollen_S),std::end(pollen_S));
    for(int value = 0; value < stigma_S.size(); value++) 
    {
		SAllele_set.insert(stigma_S[value][1]);
    }

    offspring[progenyID].slocus = stigma_S;

    // Initially set all fitness values to 1. Multiplicative fitness.
    double W_self = 1;
    double W_out = 1;
    offspring[progenyID].fitness = 1;
    offspring[progenyID].mutations = 0;
    
    // Cycle through all loci.
    for(unsigned int k = 0; k < nloci; ++k )
    {
        for(int c = 0; c < 3; ++c) // form gametes two at a time.
        {
            int count = c * 2;
            int mom = parents[gametes[count]].getGamete(wgd,k,p_rng);
            int dad = parents[gametes[count+1]].getGamete(wgd,k,p_rng);
            int offspringGenotype = mom + dad;
            offspring[progenyID].genome[k].mutations = offspringGenotype;
            
            //if(parents[gametes[count]].genome[k].mutations==1&&parents[gametes[count+1]].genome[k].mutations==2) std::cout << mom << std::endl;
            
            if(k)
            {
                if(c < 1) W_self *= offspring[progenyID].genome[k].getFitness(half_Max, ploidy);
                else if(c < 2) W_out *= offspring[progenyID].genome[k].getFitness(half_Max, ploidy);
                else {
                    offspring[progenyID].fitness *= offspring[progenyID].genome[k].getFitness(half_Max, ploidy);
                    offspring[progenyID].mutations += offspringGenotype;
                    genotypeFreqs[offspringGenotype] += 1/((double)offspring.size() * (nloci-1));
                }
            }
        }
    }
    
    // "Relative perfomance of cross types," instead of raw inbreeding depression - Ågren and Schemske 1993
    if(W_self > W_out) offspring[progenyID].ID = (W_out / W_self) - 1;
    else offspring[progenyID].ID = 1 - ( W_self / W_out );
    
    return false; // False means population did not go extinct.
}


unsigned int population::chooseParent()
{
    // The wheel of fortune model of selection. Higher fitness means a greater chance to be chosen.
    // This is a way to incorporate relative fitness into randomly choosing parents.
    
    unsigned int individualChosen = 0; //uint limit of 4,294,967,295 = max population size
    double randConverter = parents[individualChosen].fitness/fitnessSum; // Relative fitness.
    double ceiling = realDist(p_rng);
    
    while(ceiling > randConverter && individualChosen < parents.size()-1) // N-1 here because the first individual is position zero in the array
    {
        individualChosen++;
        randConverter += parents[individualChosen].fitness/fitnessSum;
    }
    
    return individualChosen;
}

void population::mutate(bool mutateS, bool reversal)
{
    if(parents.size())
    {
        unsigned int locus;
        unsigned int individual;
        
        std::poisson_distribution<> poisson( ( Mu / 2 ) * ploidy * parents.size() );
        std::uniform_int_distribution<> randomIndividual( 0, (int)parents.size() - 1 );
        std::uniform_int_distribution<> randomLocus( 0,nloci);
        std::uniform_real_distribution<> real(0,1);
        
        for(unsigned int mutations = poisson(p_rng); mutations > 0; --mutations)
        {
            locus = randomLocus(p_rng);
            individual = randomIndividual(p_rng);
            while( !mutateS && locus == 0  ) locus = randomLocus(p_rng);
            
            // We know which individual and locus gets each mutation.
            // First divide the fitness of the chosen individual by the CURRENT fitness
            // at the locus chosen. Then mutate that locus and multiply by the NEW fitness
            // at that locus.
            
            if(locus > nloci-1 ) // S-locus (locus == 1025)
            {
            	int chrom = (int)parents[individual].slocus.size() - 1;
                std::uniform_int_distribution<> randomChromosome( 0, chrom );
                if(domLevels>1)
                {
                    if( real(p_rng) > 0.5 )
                    {
                        sCounter++;
                        parents[individual].slocus[randomChromosome(p_rng)][1] = sCounter ;
                    } else
                    {
                        std::uniform_int_distribution<> dom( 0, domLevels );
                        parents[individual].slocus[randomChromosome(p_rng)][0] = dom(p_rng);
                    }
                } else
                {
                    sCounter++; 
                    parents[individual].slocus[randomChromosome(p_rng)][1] = sCounter ;
                }
            } else if(locus) // ID locus (1<=locus<=1024)
            {
                double ratio = parents[individual].genome[locus].mutations/ploidy;
                if( real(p_rng) > ratio ) // Check to see if the target allele is already a mutant. If so don't bother doing anything.
                {
                    fitnessSum -= parents[individual].fitness;
                    parents[individual].fitness /= parents[individual].genome[locus].fitness;
                    parents[individual].mutate(locus, half_Max, reversal, p_rng);
                    parents[individual].fitness *= parents[individual].genome[locus].fitness;
                    fitnessSum += parents[individual].fitness;
                }
            } else // M-locus (locus == 0)
            {
                parents[individual].mutate(locus, half_Max, reversal,p_rng);
            }
            
        }
    }
}

void population::growPop(double growthRate, double K, bool diffusion)
{
    // This function modifies the number of offspring. Take in k and r values. Calculate the pop
    // size in the next generation assuming no loss due to genetic variation then modifiy the
    // offspring array to reflect the change.
    
    double Nt = parents.size();
    
    if(Nt && Nt < K)
    {
        
        // Logistic growth function. No fractions of an offspring! Size_t truncates.
        Nt = Nt*(1 + growthRate * (1 - Nt / K));
        
        // Alternatively using lambda:
        // Nt = (lambda * Nt) \ (1 + (lambda-1) * Nt/K));
        // If lambda is used negative growth rates become possible. Right now the growthRate
        // can only drop to zero which is equivalent to a lambda of 1. If lambda is used instead
        // growth rates can be negative and mutational meltdown becomes possible. Using 
        // growthRate instead of lambda reflects a compromise between hard and soft selection.
        // To switch to lambda change every reference of growthRate back to in population.cpp
        // and main.cpp (except the log(growthRate) one)
        
        if(diffusion)
        {
            std::poisson_distribution<> poisson(Nt);
            Nt = poisson(p_rng);
        } else {
            std::bernoulli_distribution d(Nt-size_t(Nt));
            Nt = size_t(Nt) + d(p_rng);
        }
    }
    
    if( Nt >= K-1 ) Nt = K;
    
    offspring.resize( size_t(Nt), offspring[0] );
}

void population::adjustPopSize(int popSize) //
{
    // This function modifies the number of offspring possible in the next generation
    if(offspring.size()) offspring.resize( popSize, offspring[0] );
}


bool population::competePhenotypes( std::vector<std::vector<int>>& mother, std::vector<std::vector<int>>& father )
{
    // Set up two temporary storage vectors for maternal and paternal genotypes;
    std::vector<std::vector<int>> stigma_S_genotype(mother);
    std::vector<std::vector<int>> pollen_S_genotype(father);
    
    int matingRules = model;
    if(matingRules > 3) matingRules -= 3;
    
    switch ( matingRules )
    {
            
        default: // random mating
            return true;
            
        case 1: // codom model. All S-alleles have equal dominance.
            // If all alleles are co-dominante phenotype = genotype so do nothing here.
            break;
            
        case 2: // domcod model.  Pollen has dominance, pistil has codominance.
        {
            std::sort(pollen_S_genotype.begin(), pollen_S_genotype.end());
            for (unsigned int it = pollen_S_genotype.front()[0]; pollen_S_genotype.back()[0] < it; pollen_S_genotype.pop_back()) {}
            break;
        }
            
        case 3: // dom model.
        {
            std::sort(stigma_S_genotype.begin(), stigma_S_genotype.end());
            for (unsigned int it = stigma_S_genotype.front()[0]; stigma_S_genotype.back()[0] < it; stigma_S_genotype.pop_back()) {}
            
            std::sort(pollen_S_genotype.begin(), pollen_S_genotype.end());
            for (unsigned int it = pollen_S_genotype.front()[0]; pollen_S_genotype.back()[0] < it; pollen_S_genotype.pop_back()) {}
            
            break;
        }
            
            // Dom model reference: Billiard, S. Castric, V. Vekemans, X. 2007.
            // "A General Model to Explore Complex Dominance Patterns in
            // Plant Sporophytic Self-Incompatibility Systems."
    }
    
    //Compete resulting phenotypes.  Return true if reproduction is successful (no matching in phenotype), false if else.
    for (size_t n = 0; n < stigma_S_genotype.size(); n++)
    {
        for (size_t j = 0; j < pollen_S_genotype.size(); j++)
        {
            // slocus ID of 0 = disfunctional allele. First test one of the alleles for functionality.
            // If it is, test for a match. If a match is found, reproduction fails.
            if ( stigma_S_genotype[n][1] && (stigma_S_genotype[n][1] == pollen_S_genotype[j][1]) ) {
                return false;
            }
        }
    }
    //pollen_S_genotype.clear();
    return true;
}



// 2D:
// Simplify count in cross.
// Verify ploidy is needed at the individual level. If so update each offspring as they are formed in reproduce()
// Eliminate as many pointers as possible. They reduce the effectiveness of parallelization. Pointers below the population level should be fine. Make sure not to share random number engines across threads.
// update population.hpp and individual classes

// 2D2:
// FDR / SDR with S-locus. FDR / SDR with getGamete()
// Figure out why wgd freezes.

// LOOK UP FOREACH equivalent.

