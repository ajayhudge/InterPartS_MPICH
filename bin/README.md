# InterPartS Executables

This directory contains the compiled executable programs for the InterPartS simulation framework.

## Programs

### ip_prep
Preprocessor utility for InterPartS simulations. This program prepares the simulation environment and initializes necessary data structures before running the main simulation.

**Usage:**
```bash
./ip_prep
```

### ip_params2json
Utility to convert parameter files to JSON format for easier processing and validation.

**Usage:**
```bash
./ip_params2json
```

## Build Instructions

To build these executables, navigate to the `src/` directory and run:

```bash
cd ../src
make
```

## Main Simulation

After running `ip_prep`, execute the main simulation using:

```bash
mpirun -np <number_of_processors> ./iparts
```

Note: The `iparts` executable is typically located in the `src/` directory or the main build output location.
