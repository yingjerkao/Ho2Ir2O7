/* Class spin
 * Describes an Ising spin.
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

#ifndef spin_hh
#define spin_hh

#include <vec3.hh>

class tetra;

class spin {
public:
    // Possible sublattices of spins
    enum class sublattice {zero = 0, one = 1, two = 2, three = 3};

    // Directions of spin vectors
    static const vec3_int directions[4];
    
protected:
    // Sign of the Ising spin relative to standard orientation
    int m_sign;

    // Position (integer in unit-cell based coord.)
    const vec3_int m_pos;

    // Tetrahedron centres connected to this spin
    tetra* m_A;
    tetra* m_B;

    // Sublattice index
    const sublattice m_sl;
    
public:    
    //---------- CONSTRUCTION, FITTING INTO STRUCTURE -------------------------
    
    // Can only be constructed with position, direction and sublattice
    spin(const vec3_int& pos, const sublattice sl);
    
    // Virtual destructor to make sure descendants are destroyed
    virtual ~spin() {}

    // Register spin with tetrahedron centres
    void reg(tetra* t);
    // Registered tetrahedron centres
    tetra* A() const {return m_A;}
    tetra* B() const {return m_B;}    

    // Get sublattice index
    int sublat() const {return (int)m_sl;}
    
    // Get position
    vec3_int pos() const {return m_pos;}

    // Get nearest fcc lattice site
    vec3_int fcc() const {
	return (m_pos - directions[sublat()]) / 4;
    }

    //---------- BASIC SPIN OPERATIONS ---------------------------------------

    // Flip the spin
    void flip() {m_sign *= -1;}

    // Set the Ising spin
    void set(int s) {
        if (m_sign) m_sign = s;
    }

    // Get the Ising spin
    int ising() const {return m_sign;}

    // Principal direction
    vec3_int dir_plus() const {return directions[sublat()];}

    // Direction of spin
    vec3_int dir() const {return m_sign * directions[sublat()];}

    // Spin as a unit vector
    #define ISQRT3 0.577350269189625764509148780501957455647601751270126876018
    vec3 unit() const {
        return ISQRT3 * m_sign * directions[sublat()];
    }
    #undef ISQRT3
};

#endif
