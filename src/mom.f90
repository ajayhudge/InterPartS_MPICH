module mod_mom
  use mod_param
  implicit none
  private
  public momxad,momxp,momyad,momyp,momzad,momzp, momzbuoy
contains
  !
  subroutine momxad(dudt,u,v,w)
    implicit none
    real, dimension(0:,0:,0:), intent(in) :: u,v,w
    real, dimension(0:,0:,0:), intent(out) :: dudt
    integer :: im,ip,jm,jp,km,kp,i,j,k
    real :: uuip,uuim,uvjp,uvjm,uwkp,uwkm
    real :: dudxp,dudxm,dudyp,dudym,dudzp,dudzm
    !
    !$omp parallel default(none)           &
    !$omp&shared(u,v,w,dudt)               &
    !$omp&private(i,j,k,ip,jp,kp,im,jm,km) &
    !$omp&private(uuip,uuim,uvjp,uvjm,uwkp,uwkm,dudxp,dudxm,dudyp,dudym,dudzp,dudzm)
    !$omp do
    do k=1,kmax
       do j=1,jmax
          do i=1,imax
             ip = i + 1
             jp = j + 1
             kp = k + 1
             im = i - 1
             jm = j - 1
             km = k - 1
             uuip  = 0.25 * ( U(ip,j,k)+U(i,j,k) )*( U(ip,j,k)+U(i,j,k)  )
             uuim  = 0.25 * ( U(im,j,k)+U(i,j,k) )*( U(im,j,k)+U(i,j,k)  )
             uvjp  = 0.25 * ( U(i,jp,k)+U(i,j,k) )*( V(ip,j,k)+V(i,j,k)  )
             uvjm  = 0.25 * ( U(i,jm,k)+U(i,j,k) )*( V(ip,jm,k)+V(i,jm,k))
             uwkp  = 0.25 * ( U(i,j,kp)+U(i,j,k) )*( W(ip,j,k) +W(i,j,k) )
             uwkm  = 0.25 * ( U(i,j,km)+U(i,j,k) )*( W(ip,j,km)+W(i,j,km))
             dudxp = (U(ip,j,k)-U(i,j,k))*dxi
             dudxm = (U(i,j,k)-U(im,j,k))*dxi
             dudyp = (U(i,jp,k)-U(i,j,k))*dyi
             dudym = (U(i,j,k)-U(i,jm,k))*dyi
             dudzp = (U(i,j,kp)-U(i,j,k))*dzi
             dudzm = (U(i,j,k)-U(i,j,km))*dzi
             !
             ! Momentum balance
             !
             dudt(i,j,k) = dxi*( -uuip + uuim ) + visc*(dudxp-dudxm)*dxi + &
                           dyi*( -uvjp + uvjm ) + visc*(dudyp-dudym)*dyi + &
                           dzi*( -uwkp + uwkm ) + visc*(dudzp-dudzm)*dzi
          enddo
       enddo
    enddo
    !$omp end parallel
    !
    return
  end subroutine momxad
  !
  subroutine momxp(dudt,p)
    implicit none
    real, dimension(0:,0:,0:), intent(in) :: p
    real, dimension(0:,0:,0:), intent(out) :: dudt
    integer :: i,j,k
    integer :: ip,im
    !
    !$omp parallel default(none) &
    !$omp&shared(dudt,p)         &
    !$omp&private(i,j,k,ip,im)
    !$omp do
    do k=1,kmax
       do j=1,jmax
          do i=1,imax
             ip = i + 1
             im = i - 1
             dudt(i,j,k) = - dxi*( p(ip,j,k)-p(i,j,k) )
          enddo
       enddo
    enddo
    !$omp end parallel
    !
    return
  end subroutine momxp
  !
  subroutine momyad(dvdt,u,v,w)
    implicit none
    real, dimension(0:,0:,0:), intent(in) :: u,v,w
    real, dimension(0:,0:,0:), intent(out) :: dvdt
    integer :: im,ip,jm,jp,km,kp,i,j,k
    real :: uvip,uvim,vvjp,vvjm,wvkp,wvkm
    real :: dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm
    !
    !$omp parallel default(none)           &
    !$omp&shared(u,v,w,dvdt)               &
    !$omp&private(i,j,k,ip,jp,kp,im,jm,km) &
    !$omp&private(uvip,uvim,vvjp,vvjm,wvkp,wvkm,dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm)
    !$omp do
    do k=1,kmax
       do j=1,jmax
          do i=1,imax
             ip = i + 1
             jp = j + 1
             kp = k + 1
             im = i - 1
             jm = j - 1
             km = k - 1
             uvip  = 0.25 * (U(i ,j,k)+U(i ,jp,k))*( V(i,j,k)+V(ip,j,k) )
             uvim  = 0.25 * (U(im,j,k)+U(im,jp,k))*( V(i,j,k)+V(im,j,k) )
             vvjp  = 0.25 * (V(i,j,k )+V(i,jp,k) )*( V(i,j,k)+V(i,jp,k) )
             vvjm  = 0.25 * (V(i,j,k )+V(i,jm,k) )*( V(i,j,k)+V(i,jm,k) )
             wvkp  = 0.25 * (W(i,j,k )+W(i,jp,k) )*( V(i,j,kp)+V(i,j,k) )
             wvkm  = 0.25 * (W(i,j,km)+W(i,jp,km))*( V(i,j,km)+V(i,j,k) )
             dvdxp = (V(ip,j,k)-V(i,j,k))*dxi
             dvdxm = (V(i,j,k)-V(im,j,k))*dxi
             dvdyp = (V(i,jp,k)-V(i,j,k))*dyi
             dvdym = (V(i,j,k)-V(i,jm,k))*dyi
             dvdzp = (V(i,j,kp)-V(i,j,k))*dzi
             dvdzm = (V(i,j,k)-V(i,j,km))*dzi
             !
             ! Momentum balance
             !
             dvdt(i,j,k) = dxi*( -uvip + uvim ) + visc*(dvdxp-dvdxm)*dxi + &
                           dyi*( -vvjp + vvjm ) + visc*(dvdyp-dvdym)*dyi + &
                           dzi*( -wvkp + wvkm ) + visc*(dvdzp-dvdzm)*dzi
          enddo
       enddo
    enddo
    !$omp end parallel
    !
    return
  end subroutine momyad
  !
  subroutine momyp(dvdt,p)
    implicit none
    real, dimension(0:,0:,0:), intent(in) :: p
    real, dimension(0:,0:,0:), intent(out) :: dvdt
    integer :: i,j,k
    integer :: jp,jm
    !
    !$omp parallel default(none) &
    !$omp&shared(dvdt,p)         &
    !$omp&private(i,j,k,jp,jm)
    !$omp do
    do k=1,kmax
       do j=1,jmax
          do i=1,imax
             jp = j + 1
             jm = j - 1
             !
             ! Momentum balance
             !
             dvdt(i,j,k) = - dyi*( p(i,jp,k)-p(i,j,k) )
          enddo
       enddo
    enddo
    !$omp end parallel
    !
    return
  end subroutine momyp
  !
  subroutine momzad(dwdt,u,v,w)
    implicit none
    real, dimension(0:,0:,0:), intent(in) :: u,v,w
    real, dimension(0:,0:,0:), intent(out) :: dwdt
    integer :: im,ip,jm,jp,km,kp,i,j,k
    real :: uwip,uwim,vwjp,vwjm,wwkp,wwkm
    real :: dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm
    !
    !$omp parallel default(none)           &
    !$omp&shared(u,v,w,dwdt)               &
    !$omp&private(i,j,k,ip,jp,kp,im,jm,km) &
    !$omp&private(uwip,uwim,vwjp,vwjm,wwkp,wwkm,dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm)
    !$omp do
    do k=1,kmax
       do j=1,jmax
          do i=1,imax
             ip = i + 1
             jp = j + 1
             kp = k + 1
             im = i - 1
             jm = j - 1
             km = k - 1
             uwip  = 0.25 * ( W(i,j,k)+W(ip,j,k))*(U(i ,j,k)+U(i ,j,kp) )
             uwim  = 0.25 * ( W(i,j,k)+W(im,j,k))*(U(im,j,k)+U(im,j,kp) )
             vwjp  = 0.25 * ( W(i,j,k)+W(i,jp,k))*(V(i ,j,k)+V(i,j ,kp) )
             vwjm  = 0.25 * ( W(i,j,k)+W(i,jm,k))*(V(i,jm,k)+V(i,jm,kp) )
             wwkp  = 0.25 * ( W(i,j,k)+W(i,j,kp))*(W(i ,j,k)+W(i,j,kp ) )
             wwkm  = 0.25 * ( W(i,j,k)+W(i,j,km))*(W(i ,j,k)+W(i,j,km ) )
             dwdxp = (W(ip,j,k)-W(i,j,k))*dxi
             dwdxm = (W(i,j,k)-W(im,j,k))*dxi
             dwdyp = (W(i,jp,k)-W(i,j,k))*dyi
             dwdym = (W(i,j,k)-W(i,jm,k))*dyi
             dwdzp = (W(i,j,kp)-W(i,j,k))*dzi
             dwdzm = (W(i,j,k)-W(i,j,km))*dzi
             !
             ! Momentum balance
             !
             dwdt(i,j,k) = dxi*( -uwip + uwim ) + visc*(dwdxp-dwdxm)*dxi + &
                           dyi*( -vwjp + vwjm ) + visc*(dwdyp-dwdym)*dyi + &
                           dzi*( -wwkp + wwkm ) + visc*(dwdzp-dwdzm)*dzi
          enddo
       enddo
    enddo
    !$omp end parallel
    !
    return
  end subroutine momzad
  !
  subroutine momzp(dwdt,p)
    implicit none
    real, dimension(0:,0:,0:), intent(in) :: p
    real, dimension(0:,0:,0:), intent(out) :: dwdt
    integer :: kp,km
    integer :: i,j,k
    !
    !$omp parallel default(none) &
    !$omp&shared(dwdt,p)         &
    !$omp&private(i,j,k,kp,km)
    !$omp do
    do k=1,kmax
       do j=1,jmax
          do i=1,imax
             kp = k + 1
             km = k - 1 ! (SS) second order
             !
             ! Momentum balance
             !
             dwdt(i,j,k) = - dzi*( p(i,j,kp)-p(i,j,k) )
          enddo
       enddo
    enddo
    !$omp end parallel
    !
    return
  end subroutine momzp
  !
  subroutine momzbuoy(dwdt,c)
    implicit none
    real, dimension(0:,0:,0:), intent(in) :: c
    real, dimension(0:,0:,0:), intent(out) :: dwdt
    integer :: kp,km
    integer :: i,j,k
    !
    !$omp parallel default(none) &
    !$omp&shared(dwdt,p)         &
    !$omp&private(i,j,k,kp,km)
    !$omp do
    do k=1,kmax
       do j=1,jmax
          do i=1,imax
             kp = k + 1
             km = k - 1 ! (SS) second order 
             !
             ! Momentum balance
             !
             dwdt(i,j,k) =  0.5*beta*9.81*(c(i,j,k) + c(i,j,kp)) 
          enddo
       enddo
    enddo
    !$omp end parallel
    !
    return
  end subroutine momzbuoy
  !
end module mod_mom
