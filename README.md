# COSMO NBODY

Cosmo Nbody is a simulation tool that allows you to run an number of different N-body simulations, including some with cosmological expansion.
This project is created as part of the ASTRO 5470 final project for Spring 2024.

## Quick Start Guide

Clone this repository

````
git@github.com:AlexGarcia623/Nbody.git
````

Make your own simulations parameter file (see [documentation](https://github.com/AlexGarcia623/Nbody/wiki/Simulation-Parameters) for more details on parameters)

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

See our [Github Wiki](https://github.com/AlexGarcia623/Nbody/wiki) for in-depth documentation of all aspects of this project

## Issues

Feel free to file an [issue](https://github.com/AlexGarcia623/Nbody/issues) if you find a bug in the code

## All Directories

- `TestOneParams/` - contains parameter files for test1
- `TestTwoParams/` - contains parameter files for test2
- `TestThreeParams/` - contains parameter files for test3
- `img/` - contains any images used in jupyter notebooks
- `output/` - contains simulation output, is also where simulation output is written to

## All Files Explanation

- `README.md`: Provide a high-level project overview and quick start guide
- `clear_particle_data.sh`: Clear the simulation output from the `output/` directory
- `default_simulation_params.txt`: Parameters used by simulation when no file is provided
- `main.c`: Main file that runs the Nbody integration
- `makefile`: Compiler flags for C files
- `read_params.c` and `read_params.h`: Provides C and header file with functions to read in parameters and save them.
- `Test1.ipynb`: File containing test 1
- `Test2.ipynb`: File containing test 2
- `Test3.ipynb`: File containing test 3
- `read_Cosmo_Nbody.py`: Contains scripts to read simulation outputs, see [Github Wiki](https://github.com/AlexGarcia623/Nbody/wiki/Read-Cosmo-Nbody-Python-Package) 
