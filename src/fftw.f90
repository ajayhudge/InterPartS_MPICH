module mod_fftw
  implicit none
  !
  ! parameters
  !
  integer FFTW_FORWARD,FFTW_BACKWARD
  parameter (FFTW_FORWARD=-1,FFTW_BACKWARD=1)

  integer FFTW_REAL_TO_COMPLEX,FFTW_COMPLEX_TO_REAL
  parameter (FFTW_REAL_TO_COMPLEX=-1,FFTW_COMPLEX_TO_REAL=1)

  integer FFTW_ESTIMATE,FFTW_MEASURE
  parameter (FFTW_ESTIMATE=0,FFTW_MEASURE=1)

  integer FFTW_OUT_OF_PLACE,FFTW_IN_PLACE,FFTW_USE_WISDOM
  parameter (FFTW_OUT_OF_PLACE=0)
  parameter (FFTW_IN_PLACE=8,FFTW_USE_WISDOM=16)

  integer FFTW_THREADSAFE
  parameter (FFTW_THREADSAFE=128)

  integer :: plan_r2c_x, plan_c2r_x, &
       plan_r2c_y, plan_c2r_y
  !
  private
  public init_fft,fftr2c,fftc2r,clean_fft, &
       plan_r2c_x,plan_c2r_x,plan_r2c_y,plan_c2r_y
contains
  subroutine init_fft(n,plan_r2c,plan_c2r)
    implicit none
    integer, intent(in) :: n
    integer, intent(out) :: plan_r2c,plan_c2r
    real, dimension(n) :: var
    complex, dimension(n/2+1) :: varc
    !
    call dfftw_plan_dft_r2c_1d(plan_r2c,n,var,varc,FFTW_ESTIMATE)
    call dfftw_plan_dft_c2r_1d(plan_c2r,n,varc,var,FFTW_ESTIMATE)
    !
    return
  end subroutine init_fft

  subroutine fftr2c(n,varin,varout,plan_r2c)
    implicit none
    integer, intent(in) :: n,plan_r2c
    real, intent(in), dimension(n) :: varin
    complex, intent(out), dimension(n/2+1) :: varout
    !
    call dfftw_execute_dft_r2c(plan_r2c,varin,varout)
    !
    return
  end subroutine fftr2c
  !
  subroutine fftc2r(n,varin,varout,plan_c2r)
    implicit none
    integer, intent(in) :: n,plan_c2r
    complex, intent(in), dimension(n/2+1) :: varin
    real, intent(out), dimension(n) :: varout
    !
    call dfftw_execute_dft_c2r(plan_c2r,varin,varout)
    varout(:) = varout(:)/(1.*n)
    !
    return
  end subroutine fftc2r

  subroutine clean_fft(plan_r2c,plan_c2r)
    implicit none
    integer, intent(in) :: plan_r2c,plan_c2r
    !
    call dfftw_destroy_plan(plan_r2c)
    call dfftw_destroy_plan(plan_c2r)
    !
    return
  end subroutine clean_fft
  !
end module mod_fftw
