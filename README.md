# COSMO NBODY

Cosmo Nbody is a simulation tool that allows you to run an number of different N-body simulations, including some with cosmological expansion.
This project is created as part of the ASTRO 5470 final project for Spring 2024.

## Quick Start Guide

Clone this repository

````
git@github.com:AlexGarcia623/Nbody.git
````

Make your own simulations parameter file (see documentation for more details on parameters)

````
vi my_params.txt
## Enter your params in the file (see documentation for all parameters)
````

Compile the code

````
make all
````

Run the code with your parameters file

````
./main my_params.txt
````

## Documentation

See our Github Wiki for in-depth documentation of all aspects of this project

## Issues

Feel free to file an [issue](https://github.com/AlexGarcia623/Nbody/issues) if you find a bug in the code

## All Files Use

- `README.md`: Provide a high-level project overview and quick start guide
- `clear_particle_data.sh`: Clear the simulation output from the `output/` directory
- `default_simulation_params.txt`: Parameters used by simulation when no file is provided
- `main.c`: Main file that runs the Nbody integration
- `makefile`: Compiler flags for C files
- `read_params.c` and `read_params.h`: Provides C and header file with functions to read in parameters and save them.
  
