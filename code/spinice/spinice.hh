/* Class spinice
 * Spin ice simulations with single spin flip dynamics
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

#ifndef spinice_hh
#define spinice_hh

#include <vec3.hh>
#include <random>
#include <string>
#include <cstdio>

class spin;
class tetra;
class basic_stat;
class ewald_correlator;
class ewald_rs;

class spinice {
    // FCC lattice points
    static const vec3_int FCC[4];
    
    // Number of unit cells along different directions
    const int na, nb, nc;

    // Number of spins
    const size_t N;

    // Spins, Ising components stored for speed
    spin** spins;
    int* spin_signs;

    // Tetrahedra
    tetra** tetras;

    // Whether there is a dipolar interaction
    const bool dipolar;
    
private:
    // Stores Ewald summed potential terms in real space
    ewald_rs* ewald;
    const bool ext_ewald;

    // Correlators of spins, only ever need FFT'd version
    ewald_correlator* magnetic;

    // Flags to tell if divers types of correlators are to be evaluated
    const bool need_spin_corr;

    // Random number engine
    std::mt19937 random;
    // Uniform random numbers between [0,1]
    std::uniform_real_distribution<double> dist;
    
    //---------- ROUTINES TO ASSEMBLE STRUCTURE ------------------------------
    // NB basic unit of distance is a_fcc/8
    
    /* Returns spin pointer corresponding to given position,
     * provided there is a spin there */
    spin* spin_at (vec3_int pos) const;

    /* Returns tetrahedron pointer at given position,
     * provided there is a tetrahedron centre there */
    tetra* tetra_at (vec3_int pos) const;

    // Returns valid fcc lattice points inside the box
    vec3_int fcc (size_t n) const;

    // Constructs the full spin ice structure
    void init();

    //---------- PRIVATE IMPLEMENTATION OF MC --------------------------------
    // Evaluates the local "Ising field" for a given spin
    double field (size_t i) const;

    //---------- PUBLIC INTERFACE --------------------------------------------
public:
    // External magnetic field times atomic moment (viz. in units of energy)
    vec3 B;
    
    // Spin ice Hamiltonian parameters
    double J; // nearest neighbour interaction
    double D; // dipolar interaction
    double h; // intrinsic [111] field
    // Changes demag. factor and recomputes potential terms accordingly
    void set_demag(double demag);
    
    // Constructors    
    // for nearest neighbour ice
    spinice (size_t na, size_t nb, size_t nc, double J, double h, bool magnetic);
    // for dipolar ice: potential terms read from file
    spinice (size_t na, size_t nb, size_t nc, double J, double D, double h, 
             bool magnetic, std::FILE* f);
    // for dipolar ice: potential terms computed on the spot
    spinice (size_t na, size_t nb, size_t nc, double J, double D, double h,
             bool magnetic, double demag=1./3);
    // for dipolar ice: ewald_rs object passed
    spinice (size_t na, size_t nb, size_t nc, double J, double D, double h,
             bool magnetic, ewald_rs* ew);
    // Destructor
    ~spinice();

    // Seed random number generator from preset list of seeds
    void seed (int no);
    // Resets default orientations of spins
    void reset_spins();

    // MC step
    void montecarlo (double T, bool record, basic_stat* bs = NULL);

    // Number of spin sites
    size_t n_spin() const {return N;}
    /* Read-only pointer to spin
     * Arg:
     *    n: position of spin in array spins */
    const spin* get_spin(size_t n) const {return spins[n];}
    /* Read-only pointer to spin
     * Arg:
     *    pos: position of spin [a_fcc/8]
     * Throws exception if there's no spin at given point */
    const spin* get_spin(const vec3_int& pos) const {return spin_at(pos);}

    // Number of tetrahedra
    size_t n_tetra() const {return N/2;}
    /* Read-only pointer to tetrahedron
     * Arg:
     *    n: position of tetrahedron in array tetras */
    const tetra* get_tetra(size_t n) const {return tetras[n];}
    /* Read-only pointer to tetrahedron
     * Arg:
     *    pos: position of tetrahedron [a_fcc/8]
     * Throws exception if there's no tetrahedron at given point */
    const tetra* get_tetra(const vec3_int& pos) const {return tetra_at(pos);}

    /* Saves reciprocal space spin correlators
     * Filename root given as argument is extended with ".corr" */
    void save_corr(std::string s) const;

    // Saves dipole potential terms (only available for dipolar ice)
    void save_pot(std::FILE* f) const;

    // Reset correlators
    void reset_corr();
};

#endif
