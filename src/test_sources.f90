module mod_test_sources
  use mod_param
  use mod_common
  implicit none
  private
  public tests_scalar_scheme_src

contains

  subroutine tests_scalar_scheme_src(dcdt,dwdt, time)
    implicit none
    real, intent(in) :: time
    real, dimension(0:,0:,0:), intent(inout) :: dcdt
    real, dimension(0:,0:,0:), intent(inout) :: dwdt
    real :: src,artificial_buoy
    integer :: i,j,k
    real :: coorx, coory, coorz
    real :: rN, rL, theta, Edecay

    do k=1,kmax
      coorz = zc(k)/lz
      do j=1,jmax
        do i=1,imax
          ! Manufactured solution components
          coorx = xc(i)/lx
          coory = yc(j)/ly
          Edecay = exp(-8.0*pi*pi*visc*time)
          theta = sin(pi*coorx)*sin(pi*coory)*sin(pi*coorz)
          rL = -3.0*pi*pi*(visc/Pran) * theta
          rN = pi * sin(pi*yc(j))* ( &
               sin(2.0*pi*xb(i))*cos(2.0*pi*zc(k))*cos(pi*xc(i))*sin(pi*zc(k)) &
             - cos(pi*zc(k))*sin(pi*xc(i))*sin(2.0*pi*zb(k))*cos(2.0*pi*xc(i)) &
             ) * Edecay
          src = rN + rL
          artificial_buoy = theta*9.81*beta
          dcdt(i,j,k) = dcdt(i,j,k) + src
          dwdt(i,j,k) = dwdt(i,j,k) - artificial_buoy
        end do
      end do    
    end do
  end subroutine tests_scalar_scheme_src

end module mod_test_sources