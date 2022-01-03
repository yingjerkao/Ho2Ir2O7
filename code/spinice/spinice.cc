/* Implementation of class spinice
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

#include "spinice.hh"
#include "spin.hh"
#include "tetra.hh"
#include "basic_stat.hh"
#include "ewald_rs.hh"
#include "ewald_correlator.hh"
#include <vec3.hh>
#include <misc.hh>
#include <cmath>
#include <string>
#include <cstdio>
// Random numbers
#include "seed.hh"
#include <random>
#include <ctime>

using namespace std;


//-------------- ASSEMBLING STRUCTURE -----------------------------------------

// Vectors to FCC lattice sites in a cubic unit cell
const vec3_int spinice::FCC[4] = {
    vec3_int(),
    vec3_int(0,4,4),
    vec3_int(4,0,4),
    vec3_int(4,4,0)};

/* Generate fcc lattice points inside rectangular box
 * Order: pass unit cells by increasing x, y, and z in precedence
 *        within each unit cell, go in order [000]/8; [044]/8; [404]/8; [440]/8
 */
vec3_int spinice::fcc(size_t n) const {
    // Separate various indices
    int incell = n % 4;  n /= 4;
    int z      = n % nc; n /= nc;
    int y      = n % nb; n /= nb;
    int x      = n % na; // hopefully this is a trivial step

    vec3_int v(8*x,8*y,8*z);
    return v + FCC[incell];
}

/* Returns the spin sitting at some point or throws an error if there is
 * nothing there.
 * Order in array spins: fcc(n)+e[i] is at spins[4*n+i] */
spin* spinice::spin_at(vec3_int vec) const {
    int sl = -1;
    for (int i=0; i<4; i++) 
        if (alldiv(vec - spin::directions[i], 4)) {
            sl = i;
            break;
        }
    if (sl<0) throw "Invalid position of spin";
    vec -= spin::directions[sl];

    int inc = -1;
    for (int i=0; i<4; i++) 
        if (alldiv(vec-FCC[i], 8)) {
            inc = i;
            break;
        }
    if (inc<0) throw "Invalid position of fcc site";
    vec -= FCC[inc];
    
    int x = mod(vec[0]/8, na);
    int y = mod(vec[1]/8, nb);
    int z = mod(vec[2]/8, nc);
    int cell = (x*nb+y)*nc+z;

    return spins[16*cell+4*inc+sl];
}

/* Returns the spin sitting at some point or throws an error if there is
 * nothing there.
 * Order in array spins: fcc(n)+t[i] is at tetras[2*n+i] */
tetra* spinice::tetra_at(vec3_int vec) const {
    int sl = -1;
    for (int i=0; i<2; i++) 
        if (alldiv(vec-tetra::directions[i], 4)) {
            sl = i;
            break;
        }
    if (sl<0) throw "Invalid position of tetrahedron";
    vec -= tetra::directions[sl];

    int inc = -1;
    for (int i=0; i<4; i++) 
        if (alldiv(vec-FCC[i], 8)) {
            inc = i;
            break;
        }
    if (inc<0) throw "Invalid position of fcc site";
    vec -= FCC[inc];
        
    int x = mod(vec[0]/8, na);
    int y = mod(vec[1]/8, nb);
    int z = mod(vec[2]/8, nc);
    int cell = (x*nb+y)*nc+z;

    return tetras[8*cell+2*inc+sl];
}


//-------------- CONSTRUCTOR, DESTRUCTOR --------------------------------------

void spinice::init() {
    spins = new spin*[N];
    spin_signs = new int[N];
    tetras = new tetra*[N/2];

    // Creating spins and tetrahedra
    for(size_t ic = 0; ic < N/4; ic++) {
        for(int ssl = 0; ssl < 4; ssl++) {
            spins[4*ic+ssl] = new spin(fcc(ic) + spin::directions[ssl],
                                              (spin::sublattice)ssl);
            spin_signs[4*ic+ssl] = spins[4*ic+ssl] -> ising();
        }
        tetras[2*ic]   = new tetra(fcc(ic) + tetra::directions[0],
                                   tetra::sublattice::A);
        tetras[2*ic+1] = new tetra(fcc(ic) + tetra::directions[1],
                                   tetra::sublattice::B);
    }

    // Registering tetrahedra with neighbouring spins
    for(size_t i = 0; i < N/2; i++) {
        tetra* t = tetras[i];
        for(int j = 0; j < 4; j++) {
            /* t[i]->get_sl() evaluates to +/- 1,
               multiplying SL A based e[j] gives distance from actual centre */
            vec3_int r = t->pos() + t->sublat() * spin::directions[j];
            spin* s = spin_at(r);
            s->reg(t);
            t->reg(s);
        }
    }

    // Creates magnetic and electric correlator objects as necessary
    if (need_spin_corr) {
        this->magnetic = new ewald_correlator(2*na, 2*nb, 2*nc, false);
    } else {
        this->magnetic = NULL;
    }
}

spinice::spinice (size_t na, size_t nb, size_t nc, double J, double h, 
                  bool magnetic):
    na(na),
    nb(nb),
    nc(nc),
    N(16lu*na*nb*nc),
    dipolar(false),
    ext_ewald(true),
    J(J),
    D(0),
    h(h),
    need_spin_corr(magnetic)
{
    init();
    ewald = NULL;
}

spinice::spinice (size_t na, size_t nb, size_t nc, double J, double D, double h,
                  bool magnetic, std::FILE* f):
    na(na),
    nb(nb),
    nc(nc),
    N(16lu*na*nb*nc),
    dipolar(true),
    ext_ewald(false),
    J(J),
    D(D),
    h(h),
    need_spin_corr(magnetic)
{
    init();
    ewald = new ewald_rs(na, nb, nc, f);
}

spinice::spinice (size_t na, size_t nb, size_t nc, double J, double D, double h,
                  bool magnetic, double demag):
    na(na),
    nb(nb),
    nc(nc),
    N(16lu*na*nb*nc),
    dipolar(true),
    ext_ewald(false),
    J(J),
    D(D),
    h(h),
    need_spin_corr(magnetic)
{
    init();
    ewald = new ewald_rs(na, nb, nc, demag);
}

spinice::spinice (size_t na, size_t nb, size_t nc, double J, double D, double h,
                  bool magnetic, ewald_rs* ew):
    na(na),
    nb(nb),
    nc(nc),
    N(16lu*na*nb*nc),
    dipolar(true),
    ext_ewald(true),
    J(J),
    D(D),
    h(h),
    need_spin_corr(magnetic)
{
    init();
    ewald = ew;
}

spinice::~spinice() {
    for (size_t i = 0; i < N; i++) 
        delete spins[i];
    delete[] spins;
    delete[] spin_signs;

    for (size_t i = 0; i < N/2; i++)
        delete tetras[i];
    delete[] tetras;

    if (!ext_ewald)
        delete ewald;

    if (magnetic != NULL)
        delete magnetic;
}

void spinice::seed(int no) {
    random.seed(random_seed(no));
}

void spinice::set_demag(double demag) {
    ewald -> compute (demag);
}

// Reset spins into default direction
void spinice::reset_spins() {
    // sublattice 0 & 1 set to +1, sublattice 2 & 3 to -1
    for (size_t i = 0; i < N; ++i) {
        switch (spins[i]->sublat()) {
        case 0:
        case 1:
            spins[i]->set(1);
            spin_signs[i] = 1;
            break;
        default:
            spins[i]->set(-1);
            spin_signs[i] = -1;
        }
    }
}

//-------------- MONTE CARLO --------------------------------------------------

#define SQRT3 1.732050807568877293527446341505872366942805253810380628055

double spinice::field (size_t i) const {
    spin* s = spins[i];
    double retval = h + (B % s->dir_plus())/SQRT3;

    // Nearest neighbour interactions
    for (int i = 0; i < 4; ++i) {
        if (i == s->sublat()) continue;
        
        retval += J/3 * (s->A()->neighbour(i)->ising());
        retval += J/3 * (s->B()->neighbour(i)->ising());
    }

    // Dipolar interactions, if tracked
    if (dipolar) {
        int sl_i = i % 16; // index of the cubic sublattice spin i belongs to
        vec3_int cubic_i = s->fcc() / 2; // fcc indices leave no quotient
        int xi = cubic_i[0],
            yi = cubic_i[1],
            zi = cubic_i[2];
        size_t j = 0;
        for (int x = 0; x < na; ++x)
          for (int y = 0; y < nb; ++y)
            for (int z = 0; z < nc; ++z) {
                // exploit that the potential terms for the 16 consecutive spins
                // making up a cubic unit cell are contiguous in memory
                const double* pot = ewald->ptr(x-xi, y-yi, z-zi, sl_i);
                for (int jj = 0; jj < 16; ++jj)
                    retval += D * spin_signs[j+jj] * pot[jj];
                j += 16;
            }
    }

    return retval;
}

void spinice::montecarlo (double T, bool record, basic_stat* bs) {
    for (size_t i = 0; i < N; ++i) {        
        double f = field(i);

        if ( f * (spins[i]->ising()) > 0 ) {
            spins[i] -> flip();
            spin_signs[i] *= -1;
        } else if ( dist(random) < exp(2*f*(spins[i]->ising())/T) ) {
            spins[i] -> flip();
            spin_signs[i] *= -1;
        }
    }

    // Store spin correlators if needed
    if (need_spin_corr && record) {
        for (int i = 0; i < 4; ++i) {
            magnetic->clear_input();
            for (size_t n = i; n < N; n += 4) {
                vec3_int v = spins[n]->fcc();
                (*magnetic)(v[0], v[1], v[2]) = spins[n]->ising();
            }
            magnetic->set_ising(i);
        }
        magnetic->add_correlator(N);
    }

    // Record basic stats if object is provided
    if (bs) {
        bs->clear();
        for (size_t i = 0; i < N; ++i)
            bs->M += spins[i]->dir();
        for (size_t i = 0; i < N/2; ++i)
            ++(bs->Q(tetras[i]->charge(), i%2));
    }
}

//---------- SAVE TO FILE -----------------------------------------------------

static const string extension_corr(".corr");

void spinice::save_corr (string s) const {
    if (need_spin_corr) {
        FILE* f = fopen((s+extension_corr).c_str(), "wb");
        magnetic -> save_corr(f);
        fclose(f);
    }
}

void spinice::save_pot(FILE* f) const {
    if (dipolar) {
        ewald->save(f);
    } else {
        throw "Only implemented for dipolar ice";
    }
}

void spinice::reset_corr() {
    if (need_spin_corr)
        magnetic -> reset();
}
