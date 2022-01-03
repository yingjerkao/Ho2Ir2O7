/* Class tetra
 * Describes a tetrahedron centre and thus a monopole
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

#ifndef tetra_hh
#define tetra_hh

#include <vec3.hh>

class spin;

class tetra {
public:
    // Possible sublattices of tetrahedron centres
    enum class sublattice {A = 1, B = -1};

    // Relative positions of tetrahedron centres
    static const vec3_int directions[2];
    
private:
    // Neighbouring spins
    spin* neigh[4];

    // Sublattice tetrahedron belongs to
    const sublattice m_sl;
    
    // Position (integer in unit-cell based coord.)
    const vec3_int m_pos;

public:
    // Can only be constructed with position and sublattice
    tetra(const vec3_int& pos, sublattice sl): m_sl(sl), m_pos(pos) {}

    // Register with a spin
    void reg(spin* s);

    // Get sublattice index
    int sublat() const {return (int)m_sl;}
    sublattice sublat_enum() {return m_sl;}

    // Get position
    vec3_int pos() const {return m_pos;}

    // Get corresponding fcc lattice site
    vec3_int fcc() const {
	switch(m_sl) {
	case sublattice::A: return (m_pos-directions[0])/4;
	case sublattice::B: return (m_pos-directions[1])/4;
        default: throw "Invalid tetrahedron type";
	}
    }

    // Get a neighbour
    spin* neighbour(size_t n) const {return neigh[n];}

    // Monopole charge
    int charge() const;
};

#endif
