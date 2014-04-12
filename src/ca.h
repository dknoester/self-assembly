/* ca.h
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
#ifndef _CA_H_
#define _CA_H_

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <vector>

#include <ea/cmdline_interface.h>
#include <ea/datafiles/fitness.h>
#include <ea/meta_data.h>
#include <ea/line_of_descent.h>
#include <ea/mkv/markov_evolution_algorithm.h>
#include <ea/generational_models/moran_process.h>
#include <ea/selection/rank.h>
#include <ea/selection/random.h>
#include <ea/cvector.h>
#include <ea/torus.h>
#include <delay.h>
using namespace ealib;

enum objective_type { DENSITY, SYNC };

LIBEA_MD_DECL(CA_RADIUS, "self_assembly.ca.radius", int);
LIBEA_MD_DECL(CA_M, "self_assembly.ca.m", int);
LIBEA_MD_DECL(CA_N, "self_assembly.ca.n", int);
LIBEA_MD_DECL(CA_P, "self_assembly.ca.p", int);
LIBEA_MD_DECL(CA_IC_TYPE, "self_assembly.ca.initial_condition_type", int);
LIBEA_MD_DECL(CA_SAMPLES, "self_assembly.ca.samples", int);
LIBEA_MD_DECL(CA_REINFORCE, "self_assembly.ca.reinforcement", int);

#include "analysis.h"

/*! Define the EA's command-line interface.
 */
template <typename EA>
class self_assembly_cli : public cmdline_interface<EA> {
public:
    virtual void gather_options() {
        mkv::add_options(this);
        
        add_option<POPULATION_SIZE>(this);
        add_option<MORAN_REPLACEMENT_RATE_P>(this);
        add_option<RUN_UPDATES>(this);
        add_option<RUN_EPOCHS>(this);
        add_option<CHECKPOINT_PREFIX>(this);
        add_option<RNG_SEED>(this);
        add_option<RECORDING_PERIOD>(this);
        add_option<FF_INITIALIZATION_PERIOD>(this);
        add_option<DELAY_GENERATIONS>(this);
        
        add_option<CA_M>(this);
        add_option<CA_N>(this);
        add_option<CA_P>(this);
        add_option<CA_SAMPLES>(this);
        add_option<CA_RADIUS>(this);
        add_option<CA_IC_TYPE>(this);
        add_option<CA_REINFORCE>(this);
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
        add_event<lod_event>(ea);
    }
    
    //! Called before initialization (good place to calculate config options).
    virtual void before_initialization(EA& ea) {
        using namespace ealib::mkv;
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

/*! Abstract fitness function for cellular automata.
 
 Defines the initial condition matrices and consensus sequences.
 The rest (update function, state container type, etc.) is left to the subclasses.
 */
struct abstract_cellular_automata : fitness_function<unary_fitness<double>, constantS, stochasticS> {
    typedef std::vector<int> ivector_type;
    typedef boost::numeric::ublas::matrix<int> matrix_type;
    typedef boost::numeric::ublas::matrix_row<matrix_type> row_type;
    
    matrix_type _IC; //!< Matrix of initial conditions.
    ivector_type _C; //!< Consensus bit per initial condition.
    
    //! Initialize the fitness function.
    template <typename RNG, typename EA>
    void initialize(std::size_t n, RNG& rng, EA& ea) {
        using namespace std;
        
        _IC.resize(get<CA_SAMPLES>(ea), n); // initial conditions
        _C.resize(get<CA_SAMPLES>(ea)); // consensus bit
        
        switch(get<CA_IC_TYPE>(ea,0)) {
            case 0: { // uniform density + binomial
                std::size_t n=_IC.size1()/3;
                
                // 1/3 above pc:
                for(std::size_t i=0; i<n; ++i) {
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
                
                // 1/3 below pc:
                for(std::size_t i=n; i<(2*n); ++i) {
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
                
                // 1/3 at uniform probability
                for(std::size_t i=(2*n); i<_IC.size1(); ++i) {
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
};

template <typename RandomAccess>
struct reinforcement_adaptor {
    typedef typename RandomAccess::value_type value_type;
    
    reinforcement_adaptor(RandomAccess& r, std::size_t n, value_type pos, value_type neg) : _r(r), _n(n), _pos(pos), _neg(neg) {
    }

    value_type operator[](std::size_t i) {
        assert(i < (_n + 2));

        if(i < _n) {
            return _r[i];
        }
        if(i == _n) {
            return _pos;
        } else {
            return _neg;
        }
    }
    
    RandomAccess& _r;
    std::size_t _n;
    value_type _pos, _neg;
};

#endif
