module mod_gridsize
  IMPLICIT NONE
  SAVE 

  ! grid size
  integer, parameter :: itot = 40
  integer, parameter :: jtot = 40
  integer, parameter :: ktot = 80

  ! parallelisation
  integer, parameter :: ndims = 2
  integer, dimension(ndims), parameter :: dims = (/2,2/)

  ! time
  real, parameter :: t_end = 2  ! End time of simulation

  ! particles
  integer, parameter :: np = 1

  ! SHUANG: I had to define these two parameters because nl and nltot are no longer parameters.
  ! I have therefore defined a maximum number, but this is set entirely arbitrary and there 
  ! are no actual checks whether this is done appropriately. It is likely that the particle propagation
  ! does not work at the moment but we can come back to this once we are testing particle behaviour.

  integer, parameter :: nlmax = 279 !746
  integer, parameter :: nltotmax = 1915

end module mod_gridsize
