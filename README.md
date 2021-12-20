# LIF of Magnetized UCNPs

This repository contains the codes described in [1] that model the laser-induced-fluorescence (LIF) imaging process for an electron-spin polarized ultracold neutral plasma (UCNP). It contains two parts: a combined molecular-dynamics and quantum-trajectories (MDQT) simulation (**mdqt**, C++) and a rate equation model (**spectrum-model**, Matlab) that describes the population transfer of ions between states due to laser coupling and spontaneous emission.

The conventions, units, and assumptions of each code are explained in [1]. 

## mdqt

This code is an adaptation of the MDQT code from [2] and simulates the interaction of the ions in a UCNP in a magnetic field with a single imaging laser. The MD portion of the code evolves each ion's position and velocity due to inter-ion forces derived from a Yukawa one-component plasma model. The QT portion of the code evolves the ion wave functions and velocities along the LIF-laser propagation axis according to the effective Hamiltonian, which inclues the effects from an imaging laser near resonant with the <sup>2</sup>S<sub>1/2</sub>-<sup>2</sup>P<sub>1/2</sub> transition, Zeeman shifts, and spontaneous decay from the excited <sup>2</sup>P<sub>1/2</sub> manifold to the electronic ground state and the off-resonant <sup>2</sup>D<sub>3/2</sub> state.

**The simulation code requires:**
- a modern C++ compiler, such as the GNU compiler g++, please refer to [GCC, the GNU Compiler Collection](http://gcc.org) for installation and usage (using the C++17 standard or later);
- openMP parallelization, which is either built in to the compiler (Linux) or available as an external library (Mac OS);
- the [Armadillo C++ library](http://arma.sourceforge.net) for linear algebra and scientific computing.


**After installing and checking the availability of the above prerequisites, you can download, clone, and compile the source from GitHub (this example uses g++-10 and -std=c++20):**
```
$ git clone https://github.com/vrinceanu/plasma-MDQT-simulation.git
$ cd magnetized-lif-2021/mdqt
# compilation on Linux
$ g++-10 -I source/header/ -std=c++20 -fopenmp main.cpp MDQT.cpp source/md/* source/qt/* source/utils.cpp source/PlasmaSettings.cpp source/wigner-symbols/* -o main -O3 -lm -larmadillo -lstdc++fs
# run one instance of the simulation
$ ./main -p save_folder -a 0 -s ucnp.settings
```

**The input arguments to the executable (main) corresponding to each flag are:**
- *-p*: relative path to where the simulation files are saved
- *-a*: task array - this must be an integer that is zero or greater, but no larger than the number of runs minus one (see info on *runs* parameter below)
- *-s*: .settings file that controls the program

Once compiled, the MDQT code is controlled via the *.settings* file, which contains all the relevant input parameters for the plasma, LIF laser, and magnetic fields (notation consistent with [1] and units indicated in the file). 

**Description of important program options within the *.settings* file:**
- *runs*: due to the stochastic nature of particle initialization and spontaneous decay, it is sometimes useful to average over multiple runs for the same conditions. The *runs* parameter controls how many times to run for the same initial conditions. 
- *velOpt*: choose to sample particle velocities from a *maxwellian* characterized by temperature T<sub>i</sub> or *load* the velocities from pre-saved initial conditions. The *load* option is used to start simulations with equilibrated positions and velocities.
- *posOpt*: choose to initialize positions
- *thermostat*: if *true*, an Andersen thermostat is applied to the ions. During each timestep, the ions will collide with a bath of temperature T<sub>i</sub> at frequency *nu_therm*. It will be determined randomly whether each ion collides with the bath. In the event of a collision, the ion's velocity is sampled from a Maxwell-Boltzmann distribution.
- *recordForMD*: this parameter sets the frequency with which global plasma information is saved. If *true*, the sample frequency is chosen appropriately for MD and scales with the plasma frequency. If *false*, the sample frequency is chosen appropriately for the QT algorithm.
- *turnOffMD*: *true* the inter-ion forces are set to zero *false* use the Yukawa force
- *turnOffQT*: *true* turn off the evolution of ion wavefunctions *false* do not

A note on using multiple runs: the program can be configured for multiple runs by setting *runs* parameter to be greater than 1. In this case, the task array (-a) input to the executable sets the run number (run number = task array). This facilitates the use of running the executable in a for loop.

A note on using the *load* option for positions and velocities: this option is used to load from pre-saved equilibrated conditions in the equilCond/ folder. The equilCond/ folder contains the initial conditions that were loaded from to produce the results in [1]. When choosing to load, the program will find the set of conditions that most closely matches the .settings file by reading in the data folders names, which are formated as TiX_nX_GeX_NX, where X is the value for each quantity (for notation see the .settings file). Each of these pre-saved conditions was obtained with a single run of the program that initialized particle positions on a lattice and velocities from a Maxwell-Boltzmann distribution at the indicated ion temperature. These simulations were run for ~50 times the inverse ion Einstein frequency using an Andersen thermostat to ensure the desired ion temperature and spatial correlations were reached. 

**Analyzing the simulation output files:**
  - The relative save path for each run is *save_folder/runX*, where save_folder is the -p input to the executable and X is the -a (task array) integer. This folder will contain several output files containing the input parameters and ensemble-averaged plasma information such as potential and kinetic energy, velocity distributions, and the density matrix. When multiple runs are implemented, this keeps each separate run folder for the same set of conditions within the same directory, which is important for analysis.
  - There is a Matlab function for loading each output file (see [mdqt/analysis/read-files](https://github.com/vrinceanu/plasma-MDQT-simulation/tree/master/magnetized-lif-2021/mdqt/analysis/read-files)). Each file describes the structure and units of each save file.
  - There are two example analysis files for analyzing the output files. [plotIonTemps.m](https://github.com/vrinceanu/plasma-MDQT-simulation/blob/master/magnetized-lif-2021/mdqt/analysis/plotIonTemps.m) plots the run-averaged ion temperatures along each spatial dimension as a function of time. [plotStatePops.m](https://github.com/vrinceanu/plasma-MDQT-simulation/blob/master/magnetized-lif-2021/mdqt/analysis/plotStatePops.m) plots the run-averaged time evolution of the ion state populations (use to analyze sets with recordForMD=false). This file compares the MDQT simulations of the ensemble-averaged state populations to that of the rate equation (RE) and RE with Krook (REK) models from [spectrum-model](https://github.com/vrinceanu/plasma-MDQT-simulation/tree/master/magnetized-lif-2021/spectrum-model). To run either file, open it and modify the *datadir* to correspond to the appropriate save_folder (note: *datadir* should not correspond to a folder named *runX*, it should contain folders named *runX*).

**Example data sets:**
  - Two example data sets are provided in [example](https://github.com/vrinceanu/plasma-MDQT-simulation/tree/master/magnetized-lif-2021/mdqt/example), each using 10 *runs*.
  - The set labeled *dih* contains an example of using the program to simulate disorder-induced heating of a UCNP initialized with T<sub>i</sub>=1 mK and random initial positions (used recordForMD=true and turnOffQT=true). [plotIonTemps.m](https://github.com/vrinceanu/plasma-MDQT-simulation/blob/master/magnetized-lif-2021/mdqt/analysis/plotIonTemps.m) is currently set up to analyze this data set.
  - The set labeled *lif* contains a simulation of the LIF process under typical experimental conditions (used recordForMD=false and both MD and QT algorithms were implemented). [plotStatePops.m](https://github.com/vrinceanu/plasma-MDQT-simulation/blob/master/magnetized-lif-2021/mdqt/analysis/plotStatePops.m) is currently set up to analyze this data set.
  
## spectrum-model

This is a Matlab code that contains the spectrum fit model described in [1]. There are three example scripts that demonstrate how to use the models for the ion state populations, the spectrum model, and the fit model.

## References

[1] G.M. Gorman, M.K. Warrens, S.J. Bradshaw, and T.C. Killian, *Laser-Induced Fluorescence Imaging of a Spin-polarized Ultracold Neutral Plasma in a Magnetic Field*, *submitted for publication to Physical Review A*.

[2] G.M. Gorman, T.K. Langin, M. K. Warrens, D. Vrinceanu and T. C. Killian, *Combined molecular-dynamics and quantum-trajectories simulation of laser-driven, collisional systems*, Phys. Rev. A 101, 012710 (2020).
