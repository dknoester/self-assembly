/* analysis.h
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
#ifndef _ANALYSIS_H_
#define _ANALYSIS_H_

#include <boost/lexical_cast.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>

#include <ea/algorithm.h>
#include <ea/analysis.h>
#include <ea/analysis/dominant.h>
#include <ea/mkv/analysis.h>
#include "ca.h"

LIBEA_ANALYSIS_TOOL(ca_dom_1000x) {
    // get the dominant:
    typename EA::iterator i=analysis::dominant(ea);
    
    datafile df("ca_dom_1000x.dat");
    df.add_field("individual").add_field("w0").add_field("w1");
    df.write(get<IND_NAME>(*i)).write(static_cast<double>(ealib::fitness(*i,ea)));
    
    put<CA_IC_TYPE>(1,ea);
    put<CA_SAMPLES>(1000,ea);
    initialize_fitness_function(ea.fitness_function(), ea);
    
    recalculate_fitness(*i,ea);
    df.write(static_cast<double>(ealib::fitness(*i,ea))).endl();
}

LIBEA_ANALYSIS_TOOL(ca_dom_scale) {
    // get the dominant:
    typename EA::iterator i=analysis::dominant(ea);
    
    datafile df("ca_dom_scale.dat");
    df.add_field("scale").add_field("w");
    
    double m=get<CA_M>(ea);
    double n=get<CA_N>(ea);
    double p=get<CA_P>(ea);
    
    for(double s=1.0; s<10.0; s+=1.0) {
        put<CA_M>(m*s,ea);
        put<CA_N>(n*s,ea);
        put<CA_P>(p*s,ea);
        initialize_fitness_function(ea.fitness_function(), ea);
        recalculate_fitness(*i,ea);
        df.write(s).write(static_cast<double>(ealib::fitness(*i,ea))).endl();
    }
}

LIBEA_ANALYSIS_TOOL(ca_dom_ko_hidden) {
    // get the dominant:
    typename EA::iterator i=analysis::dominant(ea);
    
    datafile df("ca_dom_ko_hidden.dat");
    df.add_field("individual").add_field("w0").add_field("w1");
    df.write(get<IND_NAME>(*i)).write(static_cast<double>(ealib::fitness(*i,ea)));
    
    put<CA_IC_TYPE>(1,ea);
    put<CA_SAMPLES>(1000,ea);
    put<CA_KO_HIDDEN>(1,ea);

    initialize_fitness_function(ea.fitness_function(), ea);
    recalculate_fitness(*i,ea);
    df.write(static_cast<double>(ealib::fitness(*i,ea))).endl();
}

LIBEA_ANALYSIS_TOOL(ca_dom_adapt) {
    // get the dominant:
    typename EA::iterator i=analysis::dominant(ea);
    
    datafile df("ca_dom_adapt.dat");
    df.add_field("individual").add_field("w0").add_field("w1");
    df.write(get<IND_NAME>(*i)).write(static_cast<double>(ealib::fitness(*i,ea)));
    
    put<CA_IC_TYPE>(1,ea);
    put<CA_SAMPLES>(1000,ea);
    put<CA_REINFORCE>(0,ea);
    initialize_fitness_function(ea.fitness_function(), ea);
    
    recalculate_fitness(*i,ea);
    df.write(static_cast<double>(ealib::fitness(*i,ea))).endl();
}

LIBEA_ANALYSIS_TOOL(ca_all_100x) {
    datafile df("ca_all_1000x.dat");
    df.add_field("individual").add_field("w0").add_field("w1");
    
    put<CA_IC_TYPE>(1,ea);
    put<CA_SAMPLES>(100,ea);
    initialize_fitness_function(ea.fitness_function(), ea);
    
    for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
        df.write(get<IND_NAME>(*i)).write(static_cast<double>(ealib::fitness(*i,ea)));
        recalculate_fitness(*i,ea);
        df.write(static_cast<double>(ealib::fitness(*i,ea))).endl();
    }
}

LIBEA_ANALYSIS_TOOL(ca_dom_na_1000x) {
    // get the dominant:
    typename EA::iterator i=analysis::dominant(ea);
    
    datafile df("ca_dom_na_1000x.dat");
    df.add_field("individual").add_field("w0").add_field("w1").add_field("w_noadapt");
    df.write(get<IND_NAME>(*i)).write(static_cast<double>(ealib::fitness(*i,ea)));
    
    put<CA_IC_TYPE>(1,ea);
    put<CA_SAMPLES>(1000,ea);
    initialize_fitness_function(ea.fitness_function(), ea);
    recalculate_fitness(*i,ea);
    df.write(static_cast<double>(ealib::fitness(*i,ea)));
    
    put<CA_DISABLE_ADAPTATION>(1,ea);
    recalculate_fitness(*i,ea);
    df.write(static_cast<double>(ealib::fitness(*i,ea))).endl();
}


//! Callback used to apply noise to the state container.
template <typename FitnessFunction, typename RandomNumberGenerator>
struct noise_callback : public FitnessFunction::callback {
    noise_callback(RandomNumberGenerator& rng, double p, datafile& df) : _rng(rng), _p(p), _df(df) {
    }
    
    virtual bool new_state(typename FitnessFunction::state_container_type& s, int& c) {
        for(typename FitnessFunction::state_container_type::iterator i=s.begin(); i!=s.end(); ++i) {
            if(_rng.p(_p)) {
                (*i) ^= 0x01;
            }
        }
        return false;
    }
    
    RandomNumberGenerator& _rng;
    double _p;
    datafile& _df;
};

//! Analysis tool to apply noise to state machine inputs.
LIBEA_ANALYSIS_TOOL(ca_noise) {
    // get the dominant:
    typename EA::iterator ind=analysis::dominant(ea);
    
    put<CA_IC_TYPE>(1,ea);
    put<CA_SAMPLES>(100,ea);
    initialize_fitness_function(ea.fitness_function(), ea);
    
    datafile df("ca_noise.dat");
    df.add_field("noise_p").add_field("w");
    df.comment("individual: " + boost::lexical_cast<std::string>(get<IND_NAME>(*ind)));
    
    for(double noise=0.0; noise <= 1.0; noise += 0.02) {
        noise_callback<typename EA::fitness_function_type, typename EA::rng_type> cb(ea.rng(), noise, df);
        ea.fitness_function().reset_callback(&cb);
        recalculate_fitness(*ind,ea);
        df.write(noise).write(static_cast<double>(ealib::fitness(*ind,ea))).endl();
    }
}

LIBEA_ANALYSIS_TOOL(ca_dom_rule_density) {
    // get the dominant:
    typename EA::iterator i=analysis::dominant(ea);
    
    datafile df("ca_dom_rule_density.dat");
    df.add_field("individual").add_field("w").add_field("rho");
    df.write(get<IND_NAME>(*i)).write(static_cast<double>(ealib::fitness(*i,ea)));

    typename EA::phenotype_type& pt = ealib::phenotype(*i, ea);
    int shift = pt.nstates();
    assert(shift < (sizeof(unsigned long)*8));
    unsigned long n = 0x01 << shift;
    double rho = 0.0;
    
    for(unsigned long j=0; j<n; ++j) {
        for(int k=0; k<shift; ++k) {
            pt(k) = (j >> k) & 0x01;
        }
        pt.update();
        rho += pt.output(0);
    }
    
    df.write(rho / static_cast<double>(n)).endl();
}


LIBEA_ANALYSIS_TOOL(ca_dom_sampled_rule_density) {
    // get the dominant:
    typename EA::iterator i=analysis::dominant(ea);
    
    datafile df("ca_dom_sampled_rule_density.dat");
    df.add_field("individual").add_field("w").add_field("rho");
    df.write(get<IND_NAME>(*i)).write(static_cast<double>(ealib::fitness(*i,ea)));
    
    typename EA::phenotype_type& pt = ealib::phenotype(*i, ea);
    double rho = 0.0;
    const int trials=1000000;
    std::size_t n = pt.nstates();
    std::vector<std::size_t> bits(n);
    algorithm::iota(bits.begin(), bits.end());
    
    for(int c=0; c<trials; ++c) {
        std::random_shuffle(bits.begin(), bits.end(), ea.rng());
        int on = ea.rng()(n); // number of bits that will be on
        for(int j=0; j<on; ++j) {
            pt(bits[j]) = 1;
        }
        pt.update();
        rho += pt.output(0);
    }
    
    df.write(rho / static_cast<double>(trials)).endl();
}

//! Callback used to record a frame for a movie.
template <typename FitnessFunction>
struct movie_callback : public FitnessFunction::callback {
    movie_callback(datafile& df) : _df(df) {
    }
    
    virtual bool new_state(typename FitnessFunction::state_container_type& s, int& c) {
        _df.write_all(s.begin(), s.end()).endl();
        return false;
    }
    
    datafile& _df;
};

//! Analysis tool to generate a movie.
LIBEA_ANALYSIS_TOOL(ca_movie) {
    typename EA::iterator ind=analysis::dominant(ea);
    
    put<CA_IC_TYPE>(1,ea);
    put<CA_SAMPLES>(1,ea);
    
    datafile summary("ca_movie_summary.dat");
    summary.add_field("movie").add_field("w").add_field("c");
    summary.comment("individual: " + boost::lexical_cast<std::string>(get<IND_NAME>(*ind)));
    
    for(int i=0; i<10; ++i) {
        initialize_fitness_function(ea.fitness_function(), ea);
        
        datafile df("ca_movie_" + boost::lexical_cast<std::string>(i) + ".dat");
        df.comment("first line holds dimensions of world")
        .comment("there are always three dimensions, in (m,n,p) order.")
        .comment("if p == 0, then we're dealing with a 2d world")
        .comment("note that this is matrix notation: m=row, n=col, p=page")
        .comment("m==y axis, n==x axis, p==z axis")
        .comment("each subsequent line holds an entire world state in (page-)row-major order");
        df.write(get<CA_M>(ea,1)).write(get<CA_N>(ea)).write(get<CA_P>(ea,1)).endl();
        
        movie_callback<typename EA::fitness_function_type> cb(df);
        ea.fitness_function().reset_callback(&cb);
        recalculate_fitness(*ind,ea);
        summary.write(i).write(static_cast<double>(ealib::fitness(*ind,ea))).write(ea.fitness_function()._C[0]).endl();
    }
}

//! Callback used to record a frame for a movie and trigger adaptation.
template <typename FitnessFunction>
struct adaptive_movie_callback : public FitnessFunction::callback {
    adaptive_movie_callback(datafile& df) : _df(df), _t(0) {
    }
    
    virtual bool new_state(typename FitnessFunction::state_container_type& s, int& c) {
        _df.write_all(s.begin(), s.end()).endl();
        if(_t++==2) {// && (std::count(s.begin(), s.end(), c) > (0.8 * s.size()))) {
            c ^= 0x01;
//            _t = true;
            return true;
        }
        return false;
    }
    
    datafile& _df;
    int _t;
};

//! Analysis tool to generate a movie showing adaptation.
LIBEA_ANALYSIS_TOOL(ca_adaptive_movie) {
    typename EA::iterator ind=analysis::dominant(ea);
    
    put<CA_IC_TYPE>(1,ea);
    put<CA_SAMPLES>(1,ea);
    
    datafile summary("ca_adaptive_movie_summary.dat");
    summary.add_field("movie").add_field("w").add_field("c");
    summary.comment("individual: " + boost::lexical_cast<std::string>(get<IND_NAME>(*ind)));
    
    for(int i=0; i<10; ++i) {
        initialize_fitness_function(ea.fitness_function(), ea);
        
        datafile df("ca_adaptive_movie_" + boost::lexical_cast<std::string>(i) + ".dat");
        df.comment("first line holds dimensions of world")
        .comment("there are always three dimensions, in (m,n,p) order.")
        .comment("if p == 0, then we're dealing with a 2d world")
        .comment("note that this is matrix notation: m=row, n=col, p=page")
        .comment("m==y axis, n==x axis, p==z axis")
        .comment("each subsequent line holds an entire world state in (page-)row-major order");
        df.write(get<CA_M>(ea,1)).write(get<CA_N>(ea)).write(get<CA_P>(ea,1)).endl();
        
        adaptive_movie_callback<typename EA::fitness_function_type> cb(df);
        ea.fitness_function().reset_callback(&cb);
        recalculate_fitness(*ind,ea);
        summary.write(i).write(static_cast<double>(ealib::fitness(*ind,ea))).write(ea.fitness_function()._C[0]).endl();
    }
}

#endif
