/* Random seed supplier implementation
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

#include "seed.hh"

extern uint8_t seeds[]      asm("_binary_seed_bin_start");
extern uint8_t seeds_size[] asm("_binary_seed_bin_size");
extern uint8_t seeds_end[]  asm("_binary_seed_bin_end");

std::uint32_t random_seed(std::size_t no) {
    if (no < (std::size_t)((void *)seeds_size) / 4 ) 
        return ((uint32_t*)(seeds)) [no];
    else
        throw "There are only so many seeds";
}
        
