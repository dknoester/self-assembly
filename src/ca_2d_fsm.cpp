/* ca_1d_mkv.cpp
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
#include <boost/algorithm/string/predicate.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <vector>

#include <ea/evolutionary_algorithm.h>
#include <ea/generational_models/steady_state.h>
#include <ea/selection/rank.h>
#include <ea/selection/random.h>
#include <ea/cmdline_interface.h>
#include <ea/datafiles/fitness.h>
#include <ea/torus.h>
#include <ea/mkv/markov_evolution_algorithm.h>
using namespace ealib;

#include "ca.h"
#include "analysis.h"


/*! 1D MKV-based cellular automata.
 
 Nomenclature based on Mitchell, Crutchfield, and Das '96.
 */
struct mkv_cellular_automata : fitness_function<unary_fitness<double>, constantS, stochasticS> {
    typedef std::vector<int> ivector_type;
    typedef boost::numeric::ublas::matrix<int> matrix_type;
    typedef boost::numeric::ublas::matrix_row<matrix_type> row_type;
    typedef torus2<int> state_container_type;
    typedef offset_torus2<state_container_type> state_offset_type;
    typedef adaptor_torus2<state_offset_type> adaptor_type;
    
    matrix_type _IC; //!< Matrix of initial conditions.
    ivector_type _C; //!< Consensus bit per initial condition.
    
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
        using namespace std;
        _cb = 0;
        
        _IC.resize(get<CA_SAMPLES>(ea), get<CA_M>(ea) * get<CA_N>(ea)); // initial conditions
        _C.resize(get<CA_SAMPLES>(ea)); // consensus bit
        
        switch(get<CA_IC>(ea,0)) {
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
        
        const int m = get<CA_M>(ea);
        const int n = get<CA_N>(ea);
        const int r = get<CA_RADIUS>(ea);
        
        // build all the phenotypes, reset their rngs:
        std::vector<typename EA::phenotype_type> ca(m*n, ealib::phenotype(ind, ea));
        for(std::size_t i=0; i<ca.size(); ++i) {
            ca[i].reset(rng.seed());
        }
        
        state_container_type S_t(m, n, 0);
        state_container_type S_tplus1(m, n, 0);
        state_container_type* pt=&S_t;
        state_container_type* ptp1=&S_tplus1;
        
        // for each initial condition:
        for(std::size_t ic=0; ic<_IC.size1(); ++ic) {
            row_type row(_IC,ic);
            pt->fill(row.begin(), row.end());

            // for each update:
            for(int u=0; u<(2*m*n); ++u) {
                if(_cb != 0) {
                    _cb->new_state(*pt);
                }
                bool changed=false;
                for(int i=0; i<m; ++i) {
                    for(int j=0; j<n; ++j) {
                        int agent=n*i+j; // agent's index
                        state_offset_type offset(*pt, i-r, j-r); // offset torus to the upper left cell for this agent
                        adaptor_type adaptor(offset, 2*r+1, 2*r+1); // adapt the torus to the size of the agent's neighborhood
                        ca[agent].update(adaptor); // update the agent
                        (*ptp1)(i,j) = ca[agent].output(0); // get its output
                        changed = changed || ((*pt)(i,j) != (*ptp1)(i,j)); // and check to see if anything's changed
                    }
                }
                
                std::swap(pt, ptp1);
                if(!changed) {
                    break;
                }
            }
            
            // calculate fitness:
            switch(static_cast<objective_type>(get<CA_OBJECTIVE>(ea))) {
                case DENSITY: {
                    w += algorithm::all(S_t.begin(), S_t.end(), bind2nd(equal_to<int>(), _C[ic]));
                    break;
                }
                case SYNC: {
                    break;
                }
                default: {
                    throw bad_argument_exception("invalid ca objective");
                }
            }
        }
        return w / get<CA_SAMPLES>(ea);
    }
};


// Markov network evolutionary algorithm definition.
typedef markov_evolution_algorithm
< mkv_cellular_automata
, recombination::asexual
, generational_models::steady_state<selection::random<with_replacementS>, selection::rank>
> ea_type;

/*! Define the EA's command-line interface.
 */
template <typename EA>
class cli : public cmdline_interface<EA> {
public:
    virtual void gather_options() {
        mkv::add_options(this);
        
        add_option<POPULATION_SIZE>(this);
        add_option<STEADY_STATE_LAMBDA>(this);
        add_option<RUN_UPDATES>(this);
        add_option<RUN_EPOCHS>(this);
        add_option<CHECKPOINT_PREFIX>(this);
        add_option<RNG_SEED>(this);
        add_option<RECORDING_PERIOD>(this);
        add_option<FF_INITIALIZATION_PERIOD>(this);
        
        add_option<CA_M>(this);
        add_option<CA_N>(this);
        add_option<CA_SAMPLES>(this);
        add_option<CA_OBJECTIVE>(this);
        add_option<CA_RADIUS>(this);
        add_option<CA_IC>(this);
    }
    
    //! Define tools here.
    virtual void gather_tools() {
        add_tool<ca_dom_1000x>(this);
        add_tool<ca_all_1000x>(this);
        add_tool<ca_movie>(this);
    }

    //! Define events (e.g., datafiles) here.
    virtual void gather_events(EA& ea) {
        add_event<reinitialize_fitness_function>(ea);
        add_event<datafiles::fitness_dat>(ea);
    }
    

    //! Called before initialization (good place to calculate config options).
    virtual void before_initialization(EA& ea) {
        using namespace ealib::mkv;
        put<MKV_INPUT_N>((get<CA_RADIUS>(ea)*2+1) << 1, ea); // lshift to square the # of inputs.
        put<MKV_OUTPUT_N>(1, ea);
        
        const std::string& gates = get<MKV_GATE_TYPES>(ea);
        if(!boost::algorithm::icontains(gates, "logic")) {
            ea.config().disable(LOGIC);
        }
        if(!boost::algorithm::icontains(gates, "probabilistic")) {
            ea.config().disable(PROBABILISTIC);
        }
        if(!boost::algorithm::icontains(gates, "adaptive")) {
            ea.config().disable(ADAPTIVE);
        }
    }
};

// This macro connects the cli defined above to the main() function provided by ealib.
LIBEA_CMDLINE_INSTANCE(ea_type, cli);
