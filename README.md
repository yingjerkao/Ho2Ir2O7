## Monte Carlo simulations of pyrochlore iridates

This repository provides a "library" for Monte Carlo simulations of dipolar spin ice and related pyrochlore systems and computing thermal averages of quantities like magnetisation, emergent monopole density, and structure factors. It also includes one program file, `sweep.cc`, which computes magnetisation and monopole density in continuous field sweeps for spin ice with local [111] fields. The latter was used to generate all simulation results shown in the paper [Monopole density and antiferromagnetic domain control in spin-ice iridates](https://arxiv.org/abs/2102.04483).

### Usage

After installing the dependencies, the executable `sweep.out` can be compiled and linked using the command
```
make sweep.out
```
in the directory `code`. The simulation results presented in the paper can then be reproduced using the commands
```
./simulate.sh
./postprocess.sh
```
The latter averages the raw simulation data into 0.01 tesla bins to smoothen the curves.

Other simulation executables for other spin-ice-related projects can readily be assembled using the simulation modules available in the directory `code/spinice`

### Licence

All code in this repository is free to use and/or modify under the GNU GPL v2, see `LICENSE`. If you use the code for your academic projects, please cite the paper [Monopole density and antiferromagnetic domain control in spin-ice iridates](https://arxiv.org/abs/2102.04483).

The code depends on the [FFTW](https://www.fftw.org/) and [GSL](https://www.gnu.org/software/gsl/) libraries, also licenced under the GPL.