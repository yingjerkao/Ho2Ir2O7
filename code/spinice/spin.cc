/* Implementation of functions in class spin 
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

#include "spin.hh"
#include "tetra.hh"
#include <vec3.hh>

const vec3_int spin::directions[4] = {
    vec3_int( 1, 1, 1),
    vec3_int( 1,-1,-1),
    vec3_int(-1, 1,-1),
    vec3_int(-1,-1, 1)};

spin::spin(const vec3_int& pos, const sublattice sl):
    m_pos(pos),
    m_sl(sl)
{
    switch(m_sl) {
    case sublattice::zero:
    case sublattice::one:
        m_sign = 1;
        break;
    case sublattice::two:
    case sublattice::three:
        m_sign = -1;
        break;
    }
}

void spin::reg(tetra* t) {
    switch(t->sublat_enum()) {
    case tetra::sublattice::A: m_A = t;
        break;
    case tetra::sublattice::B: m_B = t;
        break;
    }
}
