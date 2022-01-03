# Makefile to compile and link with spinice code

CXX = g++
# CXX = icpc
CXXFLAGS = -std=c++11 -pedantic -g -O3 -Wall -no-pie
LDFLAGS = -Wl,--gc-sections

LIBS_STD = -lstdc++ -lm
LIBS_GSL = -lgsl -lgslcblas
LIBS_FFT = -lfftw3
LIBS = $(LIBS_STD) $(LIBS_GSL) $(LIBS_FFT)

SPIN_DIR = spinice
SPIN_CC = $(wildcard $(SPIN_DIR)/*.cc)
SPIN_O = $(patsubst %.cc,%.o,$(SPIN_CC)) $(SPIN_DIR)/seed_data.o

INC_PATH = -I $(SPIN_DIR) -I common

# Nontrivial header dependences
sweep.o: spinice/spinice.hh spinice/basic_stat.hh

# Generic .cc -> .o rule: note we generally don't have .hh files
%.o: %.cc
	$(CXX) $(CXXFLAGS) $(INC_PATH) $< -c -o $@

# Make an object file in another directory
$(SPIN_DIR)/%.o:
	$(MAKE) -C $(dir $@) $(notdir $@)

# Nontrivial object dependences
sweep.out: $(SPIN_O)

# Generic .o -> exec rule: uses all prerequisites (meant to be .o files)
# explicitly: 'library' object files aren't linked unless needed
%.out: %.o
	$(CXX) $(CXXFLAGS) $(INC_PATH) $(LDFLAGS) $^ $(LIBS) -o $@

# Don't allow direct compilation 
%.out: %.cc
# Don't allow removal of object files in chained implicit linking
.SECONDARY:
