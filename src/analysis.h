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

#include <ea/analysis.h>
#include <ea/analysis/dominant.h>

#include "ca.h"

LIBEA_ANALYSIS_TOOL(ca_1000x) {
    // overwrite the initial conditions embedded in the fitness function:
    put<CA_IC>(1,ea);
    put<CA_SAMPLES>(1000,ea);
    initialize_fitness_function(ea.fitness_function(), ea);
    
    // get the dominant:
    typename EA::iterator i=analysis::dominant(ea);
    
    datafile df("ca_1000x.dat");
    df.add_field("individual").add_field("w0").add_field("w1");
    df.write(get<IND_NAME>(*i)).write(static_cast<double>(fitness(*i,ea)));
    
    nullify_fitness(*i,ea);
    df.write(static_cast<double>(fitness(*i,ea))).endl();
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
    
    put<CA_IC>(1,ea);
    put<CA_SAMPLES>(10,ea);
    initialize_fitness_function(ea.fitness_function(), ea);

    for(int i=0; i<10; ++i) {
        datafile df("ca_movie_" + boost::lexical_cast<std::string>(i) + ".dat");
        df.comment("first line holds dimensions of world")
        .comment("there are always three dimensions, in (m,n,p) order.")
        .comment("if p == 0, then we're dealing with a 2d world")
        .comment("note that this is matrix notation: m=row, n=col, p=page")
        .comment("m==y axis, n==x axis, p==z axis")
        .comment("each subsequent line holds an entire world state in (page-)row-major order");
        df.write(get<CA_M>(ea)).write(get<CA_N>(ea)).write(get<CA_P>(ea,0)).endl();

        movie_callback<typename EA::fitness_function_type> cb(df);
        ea.fitness_function().reset_callback(&cb);
        recalculate_fitness(*ind,ea);
    }
}

#endif
