/* Random seed supplier implementation
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

#include "seed.hh"
#include <random>
#include <fstream>
#include <sstream>
#include <stdexcept>

// Static random number generator
static std::mt19937 rng;

void init_random(const std::string& seed_file) {
    std::ifstream file(seed_file);
    if (!file.good()) {
        throw std::runtime_error("Cannot open seed file: " + seed_file);
    }

    // Read the state from file
    std::stringstream state;
    state << file.rdbuf();
    file.close();
    
    try {
        state >> rng;
    } catch (...) {
        throw std::runtime_error("Invalid seed file format: " + seed_file);
    }
}

std::uint32_t random_seed(std::size_t no) {
    // Advance the generator based on the input number
    for (std::size_t i = 0; i < no; ++i) {
        rng.discard(1);
    }
    
    return rng();
}

bool save_random_state(const std::string& seed_file) {
    std::ofstream file(seed_file);
    if (!file.good()) {
        return false;
    }
    
    file << rng;
    file.close();
    return true;
}
        
