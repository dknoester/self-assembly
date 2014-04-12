/* ca_1d_fsm.cpp
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
#include <ea/mkv/markov_evolution_algorithm.h>
#include <ea/generational_models/moran_process.h>
#include <ea/selection/rank.h>
#include <ea/selection/random.h>
#include <ea/cvector.h>
#include <delay.h>
using namespace ealib;

#include "ca.h"


struct cellular_automata_1d : abstract_cellular_automata {
    typedef abstract_cellular_automata parent;
    typedef cvector<int> state_container_type;
    typedef cvector_offset<state_container_type> state_offset_type;

    struct callback {
        virtual void new_state(state_container_type& s) = 0;
    };
    
    callback* _cb;
    
    //! Reset the callback pointer.
    void reset_callback(callback* c=0) {
        _cb = c;
    }
    
    //! Initialize the fitness function.
    template <typename RNG, typename EA>
    void initialize(RNG& rng, EA& ea) {
        _cb = 0;
        parent::initialize(static_cast<std::size_t>(get<CA_N>(ea)), rng, ea);
    }
    
    //! Calculate fitness.
	template <typename Individual, typename RNG, typename EA>
	double operator()(Individual& ind, RNG& rng, EA& ea) {
        using namespace std;
        double w=0.0;
        
        const int n = get<CA_N>(ea);
        const int r = get<CA_RADIUS>(ea);
        
        // build all the phenotypes, reset their rngs:
        std::vector<typename EA::phenotype_type> ca(n, ealib::phenotype(ind, ea));
        for(std::size_t i=0; i<ca.size(); ++i) {
            ca[i].reset(rng.seed());
        }
        
        state_container_type S_t(n, 0);
        state_container_type S_tplus1(n, 0);
        state_container_type* pt=&S_t;
        state_container_type* ptp1=&S_tplus1;
        
        // for each initial condition:
        for(std::size_t ic=0; ic<_IC.size1(); ++ic) {
            row_type row(_IC,ic);
            std::copy(row.begin(), row.end(), pt->begin());
            
            // for each update:
            int acc=0;
            for(int u=0; u<(2*n); ++u) {
                if(_cb != 0) {
                    _cb->new_state(*pt);
                }
                bool changed=false;
                acc=0;
                for(int i=0; i<n; ++i) {
                    state_offset_type offset(*pt, i-r); // offset torus to the left cell for this agent
                    ca[i].update(offset); // update the agent
                    (*ptp1)[i] = ca[i].output(0); // get its output
                    changed = changed || ((*pt)[i] != (*ptp1)[i]); // and check to see if anything's changed
                    acc += (*ptp1)[i]; // keep an accumulation of states
                }
                std::swap(pt, ptp1);
                // early stopping:
                //  if states didn't change, or
                //  all states are 0 or 1
                if(!changed || (acc==0) || (acc==row.size())) {
                    break;
                }
            }
            
            // calculate fitness:
            w += algorithm::all(pt->begin(), pt->end(), bind2nd(equal_to<int>(), _C[ic]));
        }
        return w / get<CA_SAMPLES>(ea);
    }
};


// Markov network evolutionary algorithm definition.
typedef markov_evolution_lod_algorithm
< mean_delay<cellular_automata_1d>
, recombination::asexual
, generational_models::moran_process<selection::proportionate< >, selection::rank>
> ea_type;

// This macro connects the cli defined above to the main() function provided by ealib.
LIBEA_CMDLINE_INSTANCE(ea_type, cli);
