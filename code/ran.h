/*
 * ran_mk
 *
 * Routines for generating various kinds of random numbers.
 *
 *
 * Uses the Boost library
 * The main generator is the mt19937 generator recommended by the Boost website, which also has alternatives.
 * 
 */

#ifndef RANDOM_
#define RANDOM_


#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

extern boost::mt19937 gen;	// Declared externally in main.cpp so that it can be initialized there and called in random.cpp

typedef boost::variate_generator<boost::mt19937&, boost::uniform_int<> > rng_uniform_int_t;
typedef boost::variate_generator<boost::mt19937&, boost::uniform_real<> > rng_uniform_real_t;
typedef boost::variate_generator<boost::mt19937&, boost::exponential_distribution<> > rng_exponential_t;
typedef boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > rng_normal_t;


template<class PRNG, class Dist>
inline boost::variate_generator<PRNG&, Dist> make_gen(PRNG & rng, Dist d)
{
  return boost::variate_generator<PRNG&, Dist>(rng, d);
}

inline rng_normal_t::result_type normalrand(double mean=0.0, double var =1.0){
    
    return make_gen(gen,boost::normal_distribution<>(mean,var))();
}

inline rng_uniform_int_t::result_type randint(int lower_bound, int upper_bound)
{
//
// randint(a, b) returns a uniformly distributed random integer on [a, b]
//
    return make_gen(gen, boost::uniform_smallint<>(lower_bound, upper_bound))();
}

inline rng_uniform_real_t::result_type randreal(double lower_bound = 0, double upper_bound = 1)
{
//
// Returns a uniformly distributed random deviate.  
//	randreal() is uniform on [0, 1]
//	randreal(a, b) is uniform on [a, b].
	
// NOTE (by RG III-09): documentation on Boost says the distribution is [low, high)

    return make_gen(gen, boost::uniform_real<>(lower_bound, upper_bound))();
}

inline rng_exponential_t::result_type randexp(double lambda)
{
//
// Returns a random deviate from an exponential distribution with parameter lambda
//
    return make_gen(gen, boost::exponential_distribution<>(lambda))();
}


#endif