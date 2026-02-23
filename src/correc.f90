module mod_correc
  use mod_param
  use mod_common
  implicit none
  private
  public correc
contains
  !
  ! Corrects the velocity so that it is divergence free
  !
  subroutine correc(rkiter,p)
    real, intent(in), dimension(0:,0:,0:) :: p
    integer, intent(in) :: rkiter
    real :: factor,factori,factorj,factork
    integer :: i,j,k
    !
    factor = rkcoeffab(rkiter)*dt
    factori = factor*dxi
    factorj = factor*dyi
    factork = factor*dzi
    !
    !$omp parallel default(none) &
    !$omp&shared(factori,factorj,factork) &
    !$omp&shared(p,unew,vnew,wnew,dudt,dvdt,dwdt) &
    !$omp&private(i,j,k)
    !$omp do
    do k=1,kmax
       do j=1,jmax
          do i=1,imax
             unew(i,j,k) = dudt(i,j,k) - factori * ( p(i+1,j,k)-p(i,j,k) )
             vnew(i,j,k) = dvdt(i,j,k) - factorj * ( p(i,j+1,k)-p(i,j,k) )
             wnew(i,j,k) = dwdt(i,j,k) - factork * ( p(i,j,k+1)-p(i,j,k) )
             !      unew(i,j,k) = unew(i,j,k) - factori * ( p(i+1,j,k)-p(i,j,k) )
             !      vnew(i,j,k) = vnew(i,j,k) - factorj * ( p(i,j+1,k)-p(i,j,k) )
             !      wnew(i,j,k) = wnew(i,j,k) - factork * ( p(i,j,k+1)-p(i,j,k) )
          enddo
       enddo
    enddo
    !$omp end parallel
    !
    return
  end subroutine correc
  !
end module mod_correc
