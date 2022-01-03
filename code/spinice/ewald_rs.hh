/* Class to evaluate dipolar Ewald sums and return results in real space
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

#ifndef ewald_rs_hh
#define ewald_rs_hh

#include <vec3.hh>
#include <misc.hh>
#include <cstdio>

class ewald_rs {
    // Number of reciprocal lattice points per dimension and in total
    const int Na, Nb, Nc;
    const size_t N;
    
    /* Array of potential terms
     * potential[*][i][j] are the potential terms between a given pair of
     * sublattices __within a cubic unit cell__ (hence 16), in the order spins
     * are listed in the spinice class;
     * the first index is a row-major listing of these cubic cells
     * potential[r][i][j] is pot. between spins at R_i and r+R_j */
    double (*potential)[16][16];

    // Indexing of arrays in a row-major order as required by FFTW
    inline size_t index(int x, int y, int z) const {
        x = mod(x,Na);
        y = mod(y,Nb);
        z = mod(z,Nc);
        return (size_t(x) * Nb + y) * Nc + z;
    }
    inline size_t index(vec3_int v) const {return index(v[0],v[1],v[2]);}
    inline size_t index_fft(int x, int y, int z) const {
        x = mod(x,4*Na);
        y = mod(y,4*Nb);
        z = mod(z,4*Nc);
        return (size_t(x) * 4*Nb + y) * 4*Nc + z;
    }
    inline size_t index_fft(vec3_int v) const {
        return index_fft(v[0],v[1],v[2]);
    }

    inline double& operator()(int x, int y, int z, int i, int j) {
        return potential[index(x,y,z)][i][j];
    }

public:
    /* Constructor
     * if f is given, potential terms are read from it
     * if not, they are calculated in the constructor */
    ewald_rs(size_t Na, size_t Nb, size_t Nc, std::FILE* f);
    ewald_rs(size_t Na, size_t Nb, size_t Nc, double demag = 1.0/3.0);
    ~ewald_rs();

    // Computes potential terms for a given demag. factor
    void compute(double demag = 1.0/3.0);
    
    // Returns potential term
    inline double operator()(vec3_int v, int i, int j) const {
        return potential[index(v)][i][j];
    }

    // Returns pointer to potential term; useful for speeding up bulk access
    inline const double* ptr(int x, int y, int z, int i, int j = 0) const {
        return &(potential[index(x,y,z)][i][j]);
    }

    // Save potential terms to file (format: binary, whole array as in memory)
    void save(std::FILE* f) const;
};

#endif
