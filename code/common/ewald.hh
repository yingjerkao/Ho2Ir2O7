/* Functions related to Ewald summation.
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

#ifndef ewald_hh
#define ewald_hh

#include <gsl/gsl_sf_erf.h>

#define SQRT2 1.414213562373095048801688724209698078569671875376948073176
#define SQRTPI 1.772453850905516027298167483341145182797549456122387128213

namespace ewald {   
    /* The radius of the Gaussian blur used in Ewald summation
     * Units: cubic unit cell parameter */
    const double A = 1/SQRT2;
    const double A2 = A*A;

    // Radius (half cube side) in which real space terms are evaluated
    const int L = 10;

    /* The B and C functions appearing in the real space summation
     * the self-interaction of spins is excluded by treating r = 0 specially */
    inline double B (double r) {
        if (r == 0.0) return 0.0;
        return (erfc(r/SQRT2/A) +
                SQRT2/SQRTPI * r/A * exp(-0.5*r*r/A/A)) /
            (r*r*r);
    }
    
    inline double C (double r) {
        if (r == 0.0) return 0.0;
        return (3.0*erfc(r/SQRT2/A) +
                SQRT2/SQRTPI * r/A * (3.0 + r*r/A/A) * exp(-0.5*r*r/A/A)) /
            (r*r*r*r*r);
    }
}

#undef SQRT2
#undef SQRTPI

#endif
