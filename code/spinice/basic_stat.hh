/* Class basic_stat
 * Collects basic statistics of the spin ice sample, namely magnetisation and
 * number of monopoles of different charges and sublattices
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

#ifndef basic_stat_hh
#define basic_stat_hh

#include <vec3.hh>
#include <cstring>

struct basic_stat {
private:
    // monopole count; charge q, sublattice x stored in n_q[q+2][x]
    unsigned n_q[5][2]; 
public:
    vec3_int M; // magnetisation
    double E; // energy
    
    inline void clear() {
        M = vec3_int();
        memset(n_q, 0, 10*sizeof(unsigned));
    }

    inline unsigned& Q(int q, int x) {return n_q[q+2][x];}

    // Calculates sum of q^2 for all tetrahedra
    inline int SJ () {
        return 4*(n_q[0][0]+n_q[0][1]+n_q[4][0]+n_q[4][1]) +
            (n_q[1][0]+n_q[1][1]+n_q[3][0]+n_q[3][1]);
    }
    
    // Calculates sum of all spins (that is, 2 * sum of q on A tetrahedra)
    inline int Sh () {
        return 4*(-int(n_q[0][0])+int(n_q[4][0]))
            + 2*(-int(n_q[1][0])+int(n_q[3][0]));
    }
};

#endif
