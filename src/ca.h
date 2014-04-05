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

#include <ea/meta_data.h>

enum objective_type { DENSITY, SYNC };

LIBEA_MD_DECL(CA_RADIUS, "self_assembly.ca.radius", int);
LIBEA_MD_DECL(CA_M, "self_assembly.ca.m", int);
LIBEA_MD_DECL(CA_N, "self_assembly.ca.n", int);
LIBEA_MD_DECL(CA_SAMPLES, "self_assembly.ca.samples", int);
LIBEA_MD_DECL(CA_OBJECTIVE, "self_assembly.ca.objective", int);

#endif
