/* Implementation of class ewald_correlator.
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

#include "ewald_correlator.hh"
#include <misc.hh>
#include <complex>
#include <fftw3.h>
#include <cstdio>
// Memory operations like memset
#include <cstring>

using namespace std;

ewald_correlator::ewald_correlator (size_t Na, size_t Nb, size_t Nc, bool sum):
    n_corr(0),
    Na(Na),
    Nb(Nb),
    Nc(Nc),
    N(Na*Nb*Nc),
    sum(sum)
{
    // Creating various arrays, all of size N/2
    // Zeroing out correlator arrays
    for(size_t i = 0; i < 4; i++) {
        spin[i] = new complex<double>[N/2];
        for(size_t j = 0; j < 4; j++) {
            if (sum)
                potential[i][j] = new complex<double>[N/2];
            else
                potential[i][j] = NULL;
            correlator[i][j] = new complex<double>[N/2];
            memset(correlator[i][j], 0, N/2*sizeof(complex<double>));
        }
    }

    // Creating Fourier transform arrays using FFTW methods
    // and casting them onto complex<double>
    fft_w = fftw_alloc_complex(N);
    fft = (complex<double>*) fft_w;

    // Creating FFTW plan
    plan  = fftw_plan_dft_3d(Na, Nb, Nc, fft_w, fft_w,
                             FFTW_FORWARD, FFTW_MEASURE);
}

ewald_correlator::ewald_correlator (size_t Na, size_t Nb, size_t Nc, FILE* f):
    n_corr(0),
    Na(Na),
    Nb(Nb),
    Nc(Nc),
    N(Na*Nb*Nc),
    sum(true)
{
    // Creating various arrays, all of size N/2
    // Zeroing out correlator arrays
    // Reading potential terms from file
    for(size_t i = 0; i < 4; i++) {
        spin[i] = new complex<double>[N/2];
        for(size_t j = 0; j < 4; j++) {
            potential[i][j] = new complex<double>[N/2];
            fread(potential[i][j], sizeof(complex<double>), N/2, f);
            
            correlator[i][j] = new complex<double>[N/2];
            memset(correlator[i][j], 0, N/2*sizeof(complex<double>));
        }
    }

    // Creating Fourier transform arrays using FFTW methods
    // and casting them onto complex<double>
    fft_w = fftw_alloc_complex(N);
    fft = (complex<double>*) fft_w;

    // Creating FFTW plan
    plan  = fftw_plan_dft_3d(Na, Nb, Nc, fft_w, fft_w,
                             FFTW_FORWARD, FFTW_MEASURE);
}

ewald_correlator::~ewald_correlator() {
    // Removing various arrays
    for(size_t i = 0; i < 4; i++) {
        delete[] spin[i];
        for(size_t j = 0; j < 4; j++) {
            if (sum)
                delete[] potential[i][j];
            delete[] correlator[i][j];
        }
    }

    // Removing Fourier transform arrays
    // NB complex<double> pointers are but typecasts, no need to delete[] them
    fftw_free(fft_w);

    // Removing FFTW plan
    fftw_destroy_plan(plan);
}

// Indexing of arrays in a row-major order as required by FFTW
size_t ewald_correlator::index(int x, int y, int z) const {
    // Sanitise input, PBC
    x = mod(x, Na);
    y = mod(y, Nb);
    z = mod(z, Nc);
    return (x*Nb + y)*Nc + z;
}

void ewald_correlator::clear_input() {
    memset(fft, 0, N*sizeof(complex<double>));
}

void ewald_correlator::reset() {
    n_corr = 0; // functionality of part class reset
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            memset(correlator[i][j], 0, N/2*sizeof(complex<double>));
}

void ewald_correlator::set_potential(size_t i, size_t j) {
    if (!sum)
        throw "Cannot store potential terms.";
    // Perform Fourier transform on existing real space data
    fftw_execute(plan);
    // Copy first half to the right position
    memcpy(potential[i][j], fft, N/2*sizeof(complex<double>));
}

void ewald_correlator::set_ising(size_t i) {
    fftw_execute(plan);
    memcpy(spin[i], fft, N/2*sizeof(complex<double>));
}

/* We simply need to calculate sigma_i(q) sigma_j*(q) for all q,i,j
 * and add them to whatever already is in correlator */
void ewald_correlator::add_correlator(size_t n_sample) {
    for (size_t q=0; q<N/2; q++)
        for (size_t i=0; i<4; i++)
            for (size_t j=0; j<4; j++)
                correlator[i][j][q] += spin[i][q] * conj(spin[j][q]);
    n_corr += n_sample;
}

// Sums an array by pairwise summation
static complex<double> pairwise_sum(complex<double>* array, size_t n) {
    switch(n) {
    case 0: return 0.0;
    case 1: return *array;
    case 2: return array[0]+array[1];
    case 3: return array[0]+array[1]+array[2];
    case 4: return array[0]+array[1]+array[2]+array[3];
    default:
        size_t m = n/2;
        return pairwise_sum(array,m) + pairwise_sum(array+m,n-m);
    }
}

/* Sums V_{ij}(q) sigma_i(q) sigma_j*(q) over i,j for all q into fourier
 * Sums over q by pairwise summation, as good as FFT itself. */
double ewald_correlator::energy() {
    if (!sum)
        throw "Cannot evaluate Eawld sums";
    for (size_t q=0; q<N/2; q++) {
        fft[q] = 0.0;
        for (size_t i=0; i<4; i++)
            for (size_t j=0; j<4; j++)
                fft[q] += (potential[i][j][q] *
                           spin[i][q] * conj(spin[j][q]) );
    }
    // Sum all the numbers by pairwise summation
    return pairwise_sum(fft,N/2).real() / N;
}

// Save correlators into binary file
void ewald_correlator::save_corr(FILE* f) const {
    // array in which to divide total correlators by n_corr
    complex<double>* div = new complex<double>[N/2];
    double rec = 1.0 / n_corr;
    
    // for each array, divide each q-point by n_corr, then output
    for(size_t i=0; i<4; ++i)
        for(size_t j=0; j<4; ++j) {
            for(size_t q=0; q < N/2; ++q)
                div[q] = correlator[i][j][q] * rec;
            fwrite(div, sizeof(complex<double>), N/2, f);
        }

    // delete tmp array
    delete[] div;
}

// Save potential coeffs into binary file
void ewald_correlator::save_pot(FILE* f) const {
    for(size_t i=0; i<4; ++i)
        for(size_t j=0; j<4; ++j) 
            fwrite(potential[i][j], sizeof(complex<double>), N/2, f);
}
