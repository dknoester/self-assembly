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
    typedef cvector_offset<state_container_type> neighborhood_adaptor_type;
    typedef reinforcement_adaptor<neighborhood_adaptor_type> reinforcement_adaptor_type;
    
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
        const int nin = r*2+1;
        
        // build all the phenotypes, reset their rngs:
        std::vector<typename EA::phenotype_type> ca(n, ealib::phenotype(ind, ea));
        for(std::size_t i=0; i<ca.size(); ++i) {
            ca[i].reset(rng.seed());
            if(get<CA_DISABLE_ADAPTATION>(ea,0)) {
                ca[i].disable_adaptation();
            }
        }
        
        state_container_type S_t(n, 0);
        state_container_type S_tplus1(n, 0);
        state_container_type* pt=&S_t;
        state_container_type* ptp1=&S_tplus1;
        
        // the below are a series of adaptors that we use to collect input for each
        // cell as a subset of the entire world state.
        // we also augment this input with +/- reinforcement signals:
        neighborhood_adaptor_type neighborhood; // this lets us set the origin of a cell's neighborhood; initialized to null.
        reinforcement_adaptor_type input(neighborhood, nin); // this augments the input with +/- reinforcement signals.
        
        // for each initial condition:
        for(std::size_t ic=0; ic<_IC.size1(); ++ic) {
            for(std::size_t i=0; i<ca.size(); ++i) {
                ca[i].clear();
            }

            row_type row(_IC,ic);
            std::copy(row.begin(), row.end(), pt->begin());
            std::copy(row.begin(), row.end(), ptp1->begin());
            
            int acc=0;
            int last_acc=std::accumulate(row.begin(), row.end(), 0);
            int pos=0, neg=0;

            // for each update:
            for(int u=0; u<(2*n); ++u) {
                if(_cb != 0) {
                    _cb->new_state(*pt);
                }
                bool changed=false;
                acc=0;
                neighborhood.reset(pt); // point the neighborhood at the right state vector

                for(int i=0; i<n; ++i) {
                    int agent=i; // agent's index
                    neighborhood.reset(i-r);
                    input.reset(pos, neg);
                    ca[agent].update(input); // update the agent
                    
                    int output = ca[agent].output(0); // get its output
                    (*ptp1)[i] = output;
                    changed = changed || ((*pt)[i] != output); // and check to see if anything's changed
                    acc += output; // keep an accumulation of states
                }
                // rotate the state vector; BE SURE to use pt!
                std::swap(pt, ptp1);

                // early stopping:
                //  if states didn't change, or
                //  all states are 0 or 1
                if(!changed || (acc==0) || (acc==static_cast<int>(row.size()))) {
                    break;
                }

                // now, did we move in the right direction, compared to last time through this loop?
                pos = 0; neg = 1;
                if(_C[ic] == 1) {
                    if(acc > last_acc) {
                        pos = 1; neg = 0;
                    }
                } else {
                    if(acc < last_acc) {
                        pos = 1; neg = 0;
                    }
                }
                last_acc = acc;
            }
            // record the last state:
            if(_cb != 0) {
                _cb->new_state(*pt);
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

/*! Define the CLI for a 1D FSM.
 */
template <typename EA>
class self_assembly_1d_fsm_cli : public self_assembly_cli<EA> {
public:
    typedef self_assembly_cli<EA> parent;
    
    //! Called before initialization (good place to calculate config options).
    virtual void before_initialization(EA& ea) {
        using namespace ealib::mkv;
        int nin=get<CA_RADIUS>(ea)*2+1;
        if(get<CA_REINFORCE>(ea,0)) {
            nin += 2;
        }
        put<MKV_INPUT_N>(nin, ea);
        parent::before_initialization(ea);
    }

};

// This macro connects the cli defined above to the main() function provided by ealib.
LIBEA_CMDLINE_INSTANCE(ea_type, self_assembly_1d_fsm_cli);
