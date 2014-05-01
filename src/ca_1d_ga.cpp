/* ca_1d_ga.cpp
 *
 * This file is part of Self-Assembly.
 *
 * Copyright 2014 David B. Knoester.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <vector>

#include <ea/evolutionary_algorithm.h>
#include <ea/generational_models/steady_state.h>
#include <ea/selection/rank.h>
#include <ea/selection/random.h>
#include <ea/selection/proportionate.h>
#include <ea/representations/bitstring.h>
#include <ea/fitness_functions/all_ones.h>
#include <ea/cmdline_interface.h>
#include <ea/datafiles/fitness.h>
#include <ea/cvector.h>
using namespace ealib;

#include "ca.h"


/*! 1D GA-based cellular automata.
 
 This is a recapitulation of the work done by Mitchell, Cruthfield, and Das.
 It is not, by itself, interesting.
 */
struct ga_cellular_automata : fitness_function<unary_fitness<double>, constantS, stochasticS> {
    typedef std::vector<int> ivector_type;
    typedef boost::numeric::ublas::matrix<int> matrix_type;
    typedef boost::numeric::ublas::matrix_row<matrix_type> row_type;
    typedef cvector<int> state_container_type;
    
    matrix_type _IC; //!< Matrix of initial conditions.
    ivector_type _C; //!< Consensus bit per initial condition.
    
    struct callback {
        virtual bool new_state(state_container_type& s, int& c) = 0;
    };
    
    callback* _cb;
    
    //! Reset the callback pointer.
    void reset_callback(callback* c=0) {
        _cb = c;
    }
    
    //! Initialize the fitness function.
    template <typename RNG, typename EA>
    void initialize(RNG& rng, EA& ea) {
        using namespace std;
        _cb = 0;
        
        _IC.resize(get<CA_SAMPLES>(ea), get<CA_N>(ea)); // initial conditions
        _C.resize(get<CA_SAMPLES>(ea)); // consensus bit
        
        switch(get<CA_IC_TYPE>(ea,0)) {
            case 0: { // uniform density
                // 1/2 above pc:
                for(std::size_t i=0; i<_IC.size1()/2; ++i) {
                    row_type r(_IC,i);
                    std::size_t x=(0.5 + ea.rng().p()/2.0) * r.size();
                    for(std::size_t j=0; j<x; ++j) {
                        r[j] = 1;
                    }
                    for(std::size_t j=x; j<r.size(); ++j) {
                        r[j] = 0;
                    }
                    std::random_shuffle(r.begin(), r.end(), ea.rng());
                    _C[i] = 1;
                }
                
                // 1/2 below pc:
                for(std::size_t i=_IC.size1()/2; i<_IC.size1(); ++i) {
                    row_type r(_IC,i);
                    std::size_t x=(ea.rng().p()/2.0) * r.size();
                    for(std::size_t j=0; j<x; ++j) {
                        r[j] = 1;
                    }
                    for(std::size_t j=x; j<r.size(); ++j) {
                        r[j] = 0;
                    }
                    std::random_shuffle(r.begin(), r.end(), ea.rng());
                    _C[i] = 0;
                }
                break;
            }
            case 1: { // uniform probability
                for(std::size_t i=0; i<_IC.size1(); ++i) {
                    row_type r(_IC,i);
                    _C[i] = 0;
                    for(std::size_t j=0; j<r.size(); ++j) {
                        if(ea.rng().bit()) {
                            r[j] = 1;
                            ++_C[i];
                        } else {
                            r[j] = 0;
                        }
                    }
                    _C[i] = (_C[i] > static_cast<int>((r.size()/2)));
                }
                break;
            }
        }
    }
    
    //! Calculate fitness.
	template <typename Individual, typename RNG, typename EA>
	double operator()(Individual& ind, RNG& rng, EA& ea) {
        using namespace std;
        double w=0.0;
        
        const int n = get<CA_N>(ea);
        const int r = get<CA_RADIUS>(ea);
        
        bitstring& phi=ind.repr(); // rule table
        
        // for each initial condition:
        for(std::size_t ic=0; ic<_IC.size1(); ++ic) {
            row_type row(_IC,ic);
            state_container_type S_t(row.begin(), row.end());
            state_container_type S_tplus1(row.size(), 0);
            
            // for each update:
            for(int u=0; u<(2*n); ++u) {
                if(_cb != 0) {
                    _cb->new_state(S_t, _C[ic]);
                }
                
                bool changed=false;
                // for each element in the current state:
                for(int i=0; i<n; ++i) {
                    // build the index of the rule:
                    std::size_t x=0;
                    for(int b=(2*r), j=-r; j<=r; ++j, --b) {
                        x |= ((S_t[i+j] & 0x01) << b); // shift the low-order bit of the neighboring state b bits to the left
                    }
                    S_tplus1[i] = phi[x];
                    changed = changed || (S_t[i] != S_tplus1[i]);
                }
                std::swap(S_t, S_tplus1);
                if(!changed) {
                    break;
                }
            }
            
            // calculate fitness:
            w += algorithm::all(S_t.begin(), S_t.end(), bind2nd(equal_to<int>(), _C[ic]));
        }
        return w / static_cast<double>(get<CA_SAMPLES>(ea));
    }
};


/*! Evolutionary algorithm definition.  EAs are assembled by providing a series of
 components (representation, selection type, mutation operator, etc.) as template
 parameters.
 */
typedef evolutionary_algorithm
< individual<bitstring, ga_cellular_automata>
, ancestors::flat_bitstring
, mutation::operators::per_site<mutation::site::bitflip>
, recombination::single_point_crossover
, generational_models::steady_state<selection::random<with_replacementS>, selection::rank>
> ea_type;


/*! Define the EA's command-line interface.  Ealib provides an integrated command-line
 and configuration file parser.  This class specializes that parser for this EA.
 */
template <typename EA>
class gacli : public cmdline_interface<EA> {
public:

    //! Define the options that can be parsed.
    virtual void gather_options() {
        add_option<POPULATION_SIZE>(this);
        add_option<STEADY_STATE_LAMBDA>(this);
        add_option<MUTATION_PER_SITE_P>(this);
        add_option<RUN_UPDATES>(this);
        add_option<RUN_EPOCHS>(this);
        add_option<CHECKPOINT_PREFIX>(this);
        add_option<RNG_SEED>(this);
        add_option<RECORDING_PERIOD>(this);
        add_option<FF_INITIALIZATION_PERIOD>(this);
        
        add_option<CA_N>(this);
        add_option<CA_SAMPLES>(this);
        add_option<CA_RADIUS>(this);
        add_option<CA_IC_TYPE>(this);
    }
    
    //! Define tools here.
    virtual void gather_tools() {
        add_tool<ca_dom_1000x>(this);
        add_tool<ca_movie>(this);
    }
    
    //! Define events (e.g., datafiles) here.
    virtual void gather_events(EA& ea) {
        add_event<reinitialize_fitness_function>(ea);
        add_event<datafiles::fitness_dat>(ea);
    };
    
    //! Called before initialization (good place to calculate config options).
    virtual void before_initialization(EA& ea) {
        put<REPRESENTATION_SIZE>(0x01 << (get<CA_RADIUS>(ea)*2+1), ea);
    }
};

// This macro connects the cli defined above to the main() function provided by ealib.
LIBEA_CMDLINE_INSTANCE(ea_type, gacli);
