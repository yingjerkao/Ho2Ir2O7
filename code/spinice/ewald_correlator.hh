/* Class ewald_correlator: calculates spin correlators and Ewald summed energies
 *                         using a FT representation of the potential components
 * Life cycle:
 *     fill with real space interaction energies, Fourier transform and store
 *     input Ising spins on sublattices, Fourier transform them
 *     energy calculations can be done by
 *         E = 1/(2N) sum_{q,a,b} V_ab(q) sigma_a(q) sigma_b*(q)
 *     where N is number of FCC unit cells (i.e., 2N is what is called N in code)
 *     correlator calculations are done in terms of Fourier transformed spins
 *         (should only need Ising spin correlators by lattices)
 *
 * NB Potential terms stored in V_ab(q) correspond to real space potentials
 *    V_ab(r) between spins at R_pyro(a) and r+R_pyro(b)
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

#ifndef ewald_correlator_hh
#define ewald_correlator_hh

#include <complex>
#include <fftw3.h>

class ewald_correlator {
    // Number of stored correlators
    size_t n_corr;

    // Number of reciprocal lattice points per dimension and in total
    const size_t Na, Nb, Nc;
    const size_t N;

    // Whether or not it can actually calculate Ewald summed energies
    const bool sum;

    /* Arrays of Fourier transforms
     * Since we look at functions defined on an fcc lattice, using FFTW naively
     * creates two Brillouin zones. While this is the most efficient way of
     * doing the FFT with the library, there is no need to store but half of
     * each array, corresponding to cutting it in half along the q_x axis.
     * The stored q-points then form a (rather unconventional) reciprocal space
     * unit cell. */
    
    // Arrays to keep the potential FT
    std::complex<double>* potential[4][4];

    // Arrays to keep the spin FT
    std::complex<double>* spin[4];
    
    // Collects Ising spin correlators
    std::complex<double>* correlator[4][4];

    /* FT workhorse array
     * This must be 'full size' for the FFT to work */    
    std::complex<double>* fft;
    fftw_complex*         fft_w;

    // FFTW plan for normal Fourier transforms
    fftw_plan plan;

    // Indexing of arrays in a row-major order as required by FFTW
    size_t index(int x, int y, int z) const;

public:
    // Constructor with empty potential terms, possibly none
    ewald_correlator(size_t Na, size_t Nb, size_t Nc, bool sum);
    // Constructor reads potential terms from file
    ewald_correlator(size_t Na, size_t Nb, size_t Nc, std::FILE* f);
    // Destructor
    ~ewald_correlator();

    // Set real space input, indexing by dimension
    inline std::complex<double>& operator()(int x, int y, int z) {
        return fft[index(x,y,z)];
    }
    
    inline const
    std::complex<double>& operator()(int x, int y, int z) const {
        return fft[index(x,y,z)];
    }

    // Zero out input array
    void clear_input();

    // Zero out correlator arrays
    void reset();

    /* Set potentials corresponding to a given combination of sublattices
     * from the real space input */
    void set_potential(size_t i, size_t j);

    // Set Ising spins of a given sublattice from the real space input
    void set_ising(size_t i);

    // Get Ewald summed energy
    double energy();

    // Add current spin set-up to correlator tally
    void add_correlator(size_t n_sample);

    // Save correlators into binary file
    void save_corr(std::FILE* f) const;

    /* Save potential coefficients into binary file
     * File format is the same as that for correlators */
    void save_pot(std::FILE* f) const;
};

#endif
