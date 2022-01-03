/* Miscellaneous mathematical functions.
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

#ifndef misc_hh
#define misc_hh

#include "vec3.hh"

// Checks if all components of an integer vector is divisble by some number
inline bool alldiv(vec3_int v, int n) {
    return (!(v[0]%n)) && (!(v[1]%n)) && (!(v[2]%n));
}

// Proper modular division for potentially negative inputs
inline int mod(int a, int b) {
    if (b<0) b=-b;
    int m = a%b;
    if (m<0) m+=b;
    return m;
}

// Modular division in the range [-b/2, b/2)
inline int rem(int a, int b) {
    if (b<0) b=-b;
    int r = mod(a,b);
    if (r < (b+1)/2)
	return r;
    else
	return r-b;
}

#endif
