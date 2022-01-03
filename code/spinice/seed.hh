/* Interface to 32-bit random integers drawn from file seed.bin
 * Used as seeds to random number generators
 * seed.bin contains 16384 random bytes provided by the true random number
 * website random.org
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

#ifndef seed_hh
#define seed_hh

#include <cstdint>

std::uint32_t random_seed(std::size_t no);

#endif
