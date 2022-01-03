/* A program to generate magnetisation and monopole number history under a
 * magnetic field sweep.
 *
 * Cmdline arguments:
 *    1)   Output filename (overwritten)
 *    2)   Number of seed to be used from seed file
 *
 * The standard input has to list, in any arrangement of lines,
 *    1)    System size
 *    2)    Temperature [Kelvin]
 *    3-5)  J, D, Hloc in the Hamiltonian [Kelvin]
 *    6)    Magnetic dipole moment [Bohr magneton]
 *    7)    Demagnetisation factor
 *    8)    Magnetic field sweep rate [Oe/MC step]
 *    9)    Max. magnetic field [Oe]
 *    10-2) Direction of magnetic field (x,y,z; need not be normalised)
 *    13)   #initial equilibriation steps
 *
 * Produces a text output file: contents of each line in order:
 **   magnetic field
 **   sum of all spin vectors (multiply by mu/sqrt(3)/V for physical magn.)
 **   # of -2, -1, 0, 1, 2 charged tetrahedra on sublattice A
 **   same on sublattice B (divide by 4n^3 for density)
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

#include <vec3.hh>
#include <spinice.hh>
#include <basic_stat.hh>

#include <cstdio>
#include <cstdlib>
#include <iostream>

using namespace std;

inline void lineprint(FILE* f, double field, basic_stat& data) {
    fprintf(f, "%8.1f %5d %5d %5d", field, data.M[0], data.M[1], data.M[2]);
    for (int q = -2; q <= 2; ++q)
        fprintf(f, " %4d", data.Q(q,0));
    for (int q = -2; q <= 2; ++q)
        fprintf(f, " %4d", data.Q(q,1));
    fprintf(f, "\n");
    fflush(f);
}

int main (int argc, char** argv) {
    unsigned seed = atoi(argv[2]);
    double J, D, HLOC, MU, DEMAG, // parameters of HIO
        T, Brate, Bmax; // temperature, sweep rate, peak field
    vec3 Bdir; // direction of field
    unsigned L, burnin; // system size, initial steps
    
    // mu_B*Ga in K, field is flipped in simulations
    const double FACTOR = -6.71714e-5; 

    cin >> L >> T >> J >> D >> HLOC >> MU >> DEMAG >> Brate >> Bmax
        >> Bdir[0] >> Bdir[1] >> Bdir[2] >> burnin;

    Bdir.normalise();
    int Bstep = (int)(Bmax/Brate); // number of MC steps in one sweep

    // Set up MC simulation
    spinice s (L,L,L, J, D, HLOC, false, DEMAG);
    s.seed(seed);

    // Burn in
    s.B = vec3();
    for (unsigned i = 0; i < burnin; ++i)
        s.montecarlo(T, false);

    // Data collection and output
    basic_stat data;
    FILE* out = fopen(argv[1], "w");
    fprintf(out, "# Size: %u\n# Temperature: %f K\n# J: %f K\n# D: %f K\n",
            L,T,J,D);
    fprintf(out, "# h_loc: %f K\n# Ho moment: %f mu_B\n# Demag. factor: %f\n",
            HLOC,MU,DEMAG);
    fprintf(out, "# Burn-in samples: %u\n", burnin);
    fflush(out);

    // Sweep up 0 -> Bmax
    for (int i = 0; i < Bstep; ++i) {
        s.B = Bdir * (FACTOR*MU*Brate*i);
        s.montecarlo(T, false, &data);
        lineprint(out, Brate*i, data);
    }
    
    // Sweep down and up once
    for (int cycle = 0; cycle < 1; ++cycle) {
        for (int i = Bstep; i > -Bstep; --i) { // down
            s.B = Bdir * (FACTOR*MU*Brate*i);
            s.montecarlo(T, false, &data);
            lineprint(out, Brate*i, data);
        }
        for (int i = -Bstep; i < Bstep; ++i) { // up
            s.B = Bdir * (FACTOR*MU*Brate*i);
            s.montecarlo(T, false, &data);
            lineprint(out, Brate*i, data);
        }
    }

    fclose(out);
    return 0;
}
