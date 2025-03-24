# Makefile to compile the spinice simulator classes

CXX = g++
# CXX = icpc
CXXFLAGS = -std=c++11 -pedantic -g -O3 -Wall

INC_PATH = -I ../common -I /opt/homebrew/include

LD = ld

# Lists all object files to be made for this project
all: seed.o spin.o tetra.o ewald_correlator.o ewald_rs.o spinice.o
# Clears the directory of object files
clean:
	rm *.o
# Recompiles all object files by first removing them
remake: clean all

# Nontrivial header dependences
spin.o: tetra.hh
tetra.o: spin.hh
ewald_rs.o: spin.hh
spinice.o: spin.hh tetra.hh basic_stat.hh ewald_correlator.hh ewald_rs.hh seed.hh

# General .cc -> .o rule
%.o: %.cc %.hh 
	$(CXX) $(CXXFLAGS) $(INC_PATH) $< -c -o $@
