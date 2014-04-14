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

#include <ea/analysis.h>
#include <ea/analysis/dominant.h>

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

LIBEA_ANALYSIS_TOOL(ca_all_1000x) {
    datafile df("ca_all_1000x.dat");
    df.add_field("individual").add_field("w0").add_field("w1");
    
    put<CA_IC_TYPE>(1,ea);
    put<CA_SAMPLES>(1000,ea);
    initialize_fitness_function(ea.fitness_function(), ea);

    for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
        df.write(get<IND_NAME>(*i)).write(static_cast<double>(ealib::fitness(*i,ea)));
        recalculate_fitness(*i,ea);
        df.write(static_cast<double>(ealib::fitness(*i,ea))).endl();
    }
}

template <typename FitnessFunction>
struct movie_callback : public FitnessFunction::callback {
    movie_callback(datafile& df) : _df(df) {
    }
    
    virtual void new_state(typename FitnessFunction::state_container_type& s) {
        _df.write_all(s.begin(), s.end()).endl();
    }
    
    datafile& _df;
};

LIBEA_ANALYSIS_TOOL(ca_movie) {
    // get the dominant:
    typename EA::iterator ind=analysis::dominant(ea);
    
    put<CA_IC_TYPE>(1,ea);
    put<CA_SAMPLES>(1,ea);
    
    datafile summary("ca_movie_summary.dat");
    summary.add_field("movie").add_field("w");
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
        summary.write(i).write(static_cast<double>(ealib::fitness(*ind,ea))).endl();
    }
}

#endif
