# InterPartS Source Code

This directory contains the source code for the InterPartS simulation framework, a coupled Eulerian-Lagrangian solver for particle-laden flows.

## Directory Structure

### Main Source Files

The source directory contains the following Fortran modules (`.f90` files) and their corresponding interface files (`.i` files):

#### Core Framework
- **main.f90**: Main program entry point
- **param.f90**: Parameter definitions and global configuration
- **common.f90**: Common variables and utilities
- **initmpi.f90**: MPI initialization and setup

#### Grid and Domain
- **gridsize.f90**: Grid size computation and management
- **coordsfp.f90**: Coordinate transformation utilities (Lagrange <==> Euler)

#### Initialization
- **init.f90**: Initial conditions setup
- **initparticles.f90**: Particle initialization

#### Particle Methods
- **interp_spread.f90**: Interpolation and spreading operations (Eulerian-Lagrangian coupling)
- **correc.f90**: calculate the fluid velocity
- **mom.f90**: Momentum equations
- **scal.f90**: Scalar transport equations

#### Integration
- **rk3.f90**: Third-order Runge-Kutta time stepping
- **intgr_nwtn_eulr.f90**: Newton-Euler integration
- **intgr_over_sphere.f90**: Integration over particle surfaces

#### Physics
- **collisions.f90**: Particle collision handling
- **kernel.f90**: Kernel functions for lagrangian interpolation
- **forcing.f90**: IBM forcing
- **bound.f90**: Boundary condition handling


#### Input/Output
- **loadd.f90**: Data loading utilities
- **output.f90**: Output and result writing

#### Utilities and Diagnostics
- **chkdiv.f90**: Divergence checking
- **chkdt.f90**: Timestep checking
- **autocorr.f90**: Autocorrelation calculations
- **phase_indicator.f90**: Phase indicator calculations
- **fillps.f90**: fill the right-hand side of the Poisson equation for the correction pressure

#### Testing
- **test_sources.f90**: Test source implementations
- **tests.f90**: Test routines
- **tests_orig.f90**: Original test implementations

#### poslfp/
Lagrangian points generation code