/* Implementation of class ewald_rs
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

#include "ewald_rs.hh"
#include "spin.hh"
#include <vec3.hh>
#include <misc.hh>
#include <ewald.hh>
#include <complex>
#include <fftw3.h>
#include <cmath>
#include <cstdio>
#include <cstring>

using namespace std;
using namespace ewald;

#define SQRT2 1.414213562373095048801688724209698078569671875376948073176
#define PI 3.141592653589793238462643383279502884197169399375105820974

// Vectors to FCC lattice sites in a cubic unit cell
static const vec3_int FCC[4] = {
    vec3_int(),
    vec3_int(0,4,4),
    vec3_int(4,0,4),
    vec3_int(4,4,0)};

ewald_rs::ewald_rs(size_t na, size_t nb, size_t nc, std::FILE* f):
    Na(na), Nb(nb), Nc(nc), N(na*nb*nc)
{
    potential = new double[N][16][16];
    memset(potential, 0, 256*N*sizeof(double));
    fread(potential, 256*sizeof(double), N, f);
}

ewald_rs::ewald_rs(size_t na, size_t nb, size_t nc, double demag):
    Na(na), Nb(nb), Nc(nc), N(na*nb*nc)
{
    potential = new double[N][16][16];
    memset(potential, 0, 256*N*sizeof(double));
    compute(demag);
}

void ewald_rs::compute(double demag) {
    /* Fourier transform array and plans to evaluate reciprocal space sums
     * Real space lattice unit is a_fcc/4 */
    fftw_complex* fft_w = fftw_alloc_complex(3*64*N);
    vec3_cplx* fft = (vec3_cplx*) fft_w;
    int n[] = {4*Na,4*Nb,4*Nc};
    fftw_plan plan = fftw_plan_many_dft(3,        // rank
					n,        // dimensions
					3,        // howmany,
					fft_w, n, // in, inembed
					3, 1,     // istride, idist
					fft_w, n, // out, onembed
					3, 1,     // ostride, odist
					FFTW_FORWARD,  // sign = -1
					FFTW_ESTIMATE); // flags
    
    /* Constant prefactors of reciprocal and real space sums and surface term
     * See Ewald summation literature (and notebook for 14,19/08/19)
     * divide by 3 due to non-normalised e_i
     * Results returned in units of D for general usability */
    const double RECIP = PI / (4.0 * SQRT2) / N / 3.0;
    const double REAL = 1.0 / (16.0 * SQRT2) / 3.0;
    const double SURFACE = PI / (4.0 * SQRT2) / N / 3.0 * demag;

    // Evaluate and store `fields' due to different directions of spins
    for (int i = 0; i < 4; ++i) {
	vec3_int ei = spin::directions[i];

	// Fill out array in reciprocal space
	for (int x = -2*Na; x < 2*Na; ++x)
	    for (int y = -2*Nb; y < 2*Nb; ++y)
		for (int z = -2*Nc; z < 2*Nc; ++z) {
		    vec3 k( 2.0 * PI * x / Na,
			    2.0 * PI * y / Nb,
			    2.0 * PI * z / Nc );
		    fft[index_fft(x,y,z)] =
                        (x||y||z ?                         
                         RECIP * (ei%k)/k.len2() * exp(-0.5*A2*k.len2()) * k
                         : SURFACE * ei );
		}
	// Transform it to obtain `fields' in real space
	fftw_execute(plan);

	// Calculate reciprocal & real space term by sublattice
        // Note that each i corresponds to 4 sublattices shifted by fcc vectors
	for (int j = 0; j < 4; ++j) {
	    vec3_int ej = spin::directions[j];

            // Reciprocal space
            for (int x = 0; x < Na; ++x)
              for (int y = 0; y < Nb; ++y)
                for (int z = 0; z < Nc; ++z)
                  for (int fi = 0; fi < 4; ++fi)
                    for (int fj = 0; fj < 4; ++fj) {
                        // Position in FFT array
                        vec3_int rf = (8*vec3_int(x,y,z) + (FCC[fj] + ej) -
                                       (FCC[fi] + ei)) / 2;
                        (*this)(x,y,z, 4*fi+i, 4*fj+j) =
                            fft[index_fft(rf)].real() % ej;
                    }

            // Real space - (2L)^3 cubic unit cells
	    for (int x = -L; x < L; ++x)
              for (int y = -L; y < L; ++y)
		for (int z = -L; z < L; ++z)
                  for (int fi = 0; fi < 4; ++fi)
                    for (int fj = 0; fj < 4; ++fj) {
                        // Physical separation vector
                        vec3 r = (8*vec3_int(x,y,z) + (FCC[fj] + ej) -
                                  (FCC[fi] + ei)) / 8.0;
                        double l = r.len();

                        if (l < 1e-3) {
                            // Self-interaction is set to zero
                            (*this)(x,y,z, 4*fi+i, 4*fj+j) = 0.0;
                        } else {
                            // Otherwise add real space term
                            (*this)(x,y,z, 4*fi+i, 4*fj+j) +=
                                REAL * ((ei%ej) * B(l) - (ei%r)*(ej%r) * C(l));
                        }
                    }
        }
    }

    // Ensure there is no self-interaction
    for (int i = 0; i < 16; ++i)
        (*this)(0,0,0,i,i) = 0.0;

    fftw_free(fft_w);
    fftw_destroy_plan(plan);
}

ewald_rs::~ewald_rs() {
    delete[] potential;
}

void ewald_rs::save(std::FILE* f) const {
    fwrite(potential, 256*sizeof(double), N, f);
}
    
