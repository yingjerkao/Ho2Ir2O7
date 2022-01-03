/* Implementation of functions in class tetra 
 *
 * Copyright (C) 2022 Attila Szabó <attila.szabo@physics.ox.ac.uk>
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License,
 * version 2, as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 */

#include "tetra.hh"
#include "spin.hh"
#include <vec3.hh>

const vec3_int tetra::directions[2] = {
    vec3_int(),
    vec3_int(2,2,2)};

void tetra::reg(spin* s) {neigh[s->sublat()] = s;}

int tetra::charge() const {
    int q = (neigh[0] -> ising() +
             neigh[1] -> ising() +
             neigh[2] -> ising() +
             neigh[3] -> ising() );
    /* Ising spin is +1 if it points out of an A-tetra, thus into a B-tetra
     * => there is a sign difference between the two cases
     * Divide by 2 to get +/-1 for single monopole*/
    return q * sublat() /2;
}