/* Interface to random seed generation for random number generators
 *
 * Copyright (C) 2024 Ying-Jer Kao <yjkao@phys.ntu.edu.tw>
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
#include <cstddef>
#include <string>
#include <stdexcept>

// Default seed file name
const std::string DEFAULT_SEED_FILE = "random_seed.txt";

// Initialize the random number generator from a seed file
// Throws std::runtime_error if the file doesn't exist or is invalid
void init_random(const std::string& seed_file = DEFAULT_SEED_FILE);

// Returns a random seed based on the input number
std::uint32_t random_seed(std::size_t no);

// Save current state of the random number generator
// Returns true if save was successful, false otherwise
bool save_random_state(const std::string& seed_file = DEFAULT_SEED_FILE);

#endif
