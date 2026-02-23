module mod_param
  use decomp_2d
  use mod_gridsize
  implicit none

  
  ! DEFINTIONS 

  ! run mode 
  integer, parameter :: RUNMODE_NORMAL = 0
  integer, parameter :: RUNMODE_TEST_RK_STAGE      = 900
  integer, parameter :: RUNMODE_TEST_SCALAR_BCS    = 1000
  integer, parameter :: RUNMODE_TEST_SCALAR_SCHEME = 1001
  integer, parameter :: RUNMODE_TEST_SCALAR_ACTIVE = 1002
  integer, parameter :: RUNMODE_TEST_PAR_FREE_FALL_NO_SCALAR = 1003

  ! VARIABLES
  
  ! Domain size
  real :: lx = 1.,ly = 1.,lz = 1.
  
  ! Type of initial velocity field (see init.f90)
  ! iniu = 'cou' --> plane Couette flow
  !      = 'poi' --> plane Poiseuille flow
  !      = 'zer' --> zero velocity everywhere
  !      = 'log' --> logarithmic profile + random noise
  character(len=9) :: iniu = 'zer' !

  ! Type of initial scalar concentration field (see init.f90)
  character(len=9) :: inic = 'ZERC' ! Zero concentration everywhere

  ! Boundary conditions

  ! Velocity BCs
  character(len=15) :: bc_uvw_top_type = 'BC_NOSLIP'
  character(len=15) :: bc_uvw_bot_type = 'BC_NOSLIP'

  ! Scalar BCs
  character(len=15) :: bc_c_top_type = 'BC_DIRICHLET'
  character(len=15) :: bc_c_bot_type = 'BC_DIRICHLET'
  real :: bc_c_top_val = 0. ! Dirichlet: top boundary scalar value
  real :: bc_c_bot_val = 0. ! Dirichlet: bottom boundary scalar value

  ! Other physical variables
  real :: visc = 3.21e-3 ! Kinematic viscosity
  real :: Pran = 1.0 ! Prandtl number
  real :: beta = 0. ! Expansion coefficient of scalar
  real :: radius = 0.5 ! particle radius
  real :: ratiorho = 1.5 ! particle to fluid density ratio

  ! Particle rotation control
  logical :: particle_rotation = .false. ! set to .true. to enable rotation

  ! Run mode
  integer :: runmode = RUNMODE_TEST_SCALAR_ACTIVE

  ! TO BE ADDED TO NAMELIST (PROBABLY)

  !Output parameters
  integer :: ioutchk = 5
  integer :: iout1d = 5
  integer :: iout2d = 3000
  integer :: ioutfld = 3000

  ! DEFINE NAMELIST
  NAMELIST/parameters/ &
    lx, ly, lz, &
    iniu, &
    inic, &
    bc_uvw_top_type, &
    bc_uvw_bot_type, &
    bc_c_top_type, &
    bc_c_bot_type, &
    bc_c_top_val, &
    bc_c_bot_val, &
    visc, &
    Pran, &
    beta, &
    radius, &
    ratiorho, &
    particle_rotation, &
    ioutchk, &
    iout1d, & 
    iout2d, &
    ioutfld, &
    runmode
  ! --------------------------------------------------------------------------
  ! INTERNAL PARAMETERS

  ! Grid parameters
  integer, parameter :: it1 = itot+1, jt1 = jtot+1, kt1 = ktot+1
  integer, parameter :: imax = itot/dims(1), jmax = jtot/dims(2),kmax = ktot
  integer, parameter :: i1 = imax+1, j1 = jmax+1, k1 = kmax+1

  real, parameter :: pi = acos(-1.)
  real, parameter :: epsilon = 0.01 ! Amplitude of the initial RBC perturbation (SS)
  real, parameter :: h_inter = 0.033 ! 

  integer, parameter :: npmax = nint(min(1.*np,max(1.,5.*np/(1.*dims(1)*dims(2)))))
  integer, parameter :: nqmax = min(np+2,15) ! np+2 walls
  
  real, parameter :: gacc = -9.81 !diam/vscale**2.
  real, parameter :: gaccx = gacc*cos(pi/2.), gaccy = gacc*cos(pi/2.), &
       gaccz = gacc*sin(pi/2.)

  ! Particle parameters
  integer, parameter :: nfriendsmax = 80

  ! Timestep parameters
  real, parameter, dimension(3,2) :: rkcoeff = reshape((/ 32./60., 25./60., 45./60., 0., -17./60., -25./60. /), shape(rkcoeff))
  ! rkcoeff = [32/60 0; 25/60 -17/60; 45/60 -25/60]
  real, parameter, dimension(3) :: rkcoeffab = rkcoeff(:,1)+rkcoeff(:,2)
  ! rkcoeffab = [32/60; 8/60; 20/60]
  ! Collision parameters
  real, parameter :: Nstretch = 8, dt_estim = 0.05!0.003
  real, parameter :: en = 0.97, et = 0.10, muc = 0.

  ! Lubrication model
  real, parameter :: eps_ini_pp = 0.025, eps_sat_pp = 0.001, eps_cut_pp = 0. ! sphere/sphere
  real, parameter :: eps_ini_pw = 0.075, eps_sat_pw = 0.001, eps_cut_pw = 0. ! sphere/wall

  ! --------------------------------------------------------------------------
  ! DERIVED PARAMETERS

  ! grid
  real :: dxi, dyi, dzi
  real :: dx, dy, dz
  real, dimension(0:i1) :: xc, xb    ! local coordinates ! global coordinates (SS)
  real, dimension(0:j1) :: yc, yb
  real, dimension(0:k1) :: zc, zb

  ! mpi
  integer :: send_real, send_int ! amount of data to be send from mstr to slve: see common file

  ! Thermal diffusivity
  real :: kappa
  
  ! particles
  real :: offset
  real :: volp
  real :: mominert
  real :: lref, uref, tref
  real :: retrac

  integer :: nl, NL2, NL3, NL4, NLtot
  real :: solidity

  ! sphere/sphere
  real :: colthr_pp
  real :: meffn_ss, mefft_ss, kn_ss, kt_ss, etan_ss, etat_ss, muc_ss

  ! sphere/wall
  real :: colthr_pw
  real :: meffn_sw, mefft_sw, kn_sw, kt_sw, etan_sw, etat_sw, muc_sw, psi_crit_sw, psi_crit_ss
  
  !
  ! set parameters for the lubrication model
  !
  real :: coeff_f, coeff_t

  ! Correction with Stokes amplification factor is added for values of epsilon smaller than
  ! eps_ini_pp
  real :: a11_ini_pp, a11a_ini_pp, a11b_ini_pp, a22_ini_pp, a22a_ini_pp, a22b_ini_pp, &
          a33a_ini_pp, a33b_ini_pp, b23_ini_pp, b23a_ini_pp, b23b_ini_pp, b32a_ini_pp, &
          b32b_ini_pp, c23a_ini_pp, c23b_ini_pp, c32a_ini_pp, c32b_ini_pp, d11_ini_pp, &
          d11a_ini_pp, d11b_ini_pp, d22a_ini_pp, d22b_ini_pp, d33a_ini_pp, d33b_ini_pp
  
  ! Stokes amplification factor saturated for values of epsilon smaller than eps_sat_pp
  real :: a11_sat_pp, a11a_sat_pp, a11b_sat_pp, a22_sat_pp, a22a_sat_pp, a22b_sat_pp, & 
          a33a_sat_pp, a33b_sat_pp, b23_sat_pp, b23a_sat_pp, b23b_sat_pp, b32a_sat_pp, &
          b32b_sat_pp, c23a_sat_pp, c23b_sat_pp, c32a_sat_pp, c32b_sat_pp, d11_sat_pp, &
          d11a_sat_pp, d11b_sat_pp, d22a_sat_pp, d22b_sat_pp, d33a_sat_pp, d33b_sat_pp

  ! Correction with Stokes amplification factor is added for values of epsilon smaller than
  ! eps_ini_pw
  real :: a11_ini_pw, a22_ini_pw, a33_ini_pw, b23_ini_pw, b32_ini_pw, c23_ini_pw, &
          c32_ini_pw, d11_ini_pw, d22_ini_pw, d33_ini_pw

  ! Stokes amplification factor saturated for values of epsilon smaller than eps_sat_pw
  real :: a11_sat_pw, a22_sat_pw, a33_sat_pw, b23_sat_pw, b32_sat_pw, &
       c23_sat_pw, c32_sat_pw, d11_sat_pw, d22_sat_pw, d33_sat_pw
  !
end module mod_param
