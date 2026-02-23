module mod_scal
  use mod_param
  use mod_common
  use mod_bound
  use mpi
  implicit none
  private
  public scalad
contains
  !
  subroutine scalad(dcdt,u,v,w,c) 
   implicit none
   real, dimension(0:,0:,0:), intent(in) :: u,v,w,c
   real, dimension(0:,0:,0:), intent(out) :: dcdt
   integer :: im,ip,jm,jp,km,kp,i,j,k
   real :: ucim,ucip,vcjm,vcjp,wckm,wckp
   real :: dcdxp,dcdxm,dcdyp,dcdym,dcdzp,dcdzm
   !
   !$omp parallel default(none)           &
   !$omp&shared(u,v,w,c,dudt,dcdt)               &
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
             
             ucip  = 0.5*U(i,j,k)*( C(ip,j,k)+C(i,j,k)  )
             ucim  = 0.5*U(im,j,k)*( C(i,j,k)+C(im,j,k)  )
             vcjp  = 0.5*V(i,j,k)*( C(i,jp,k)+C(i,j,k)  )
             vcjm  = 0.5*V(i,jm,k)*( C(i,j,k)+C(i,jm,k)  )
             wckp  = 0.5*W(i,j,k)*( C(i,j,kp)+C(i,j,k)  )
             wckm  = 0.5*W(i,j,km)*( C(i,j,k)+C(i,j,km)  )

             dcdxp = (C(ip,j,k)-C(i,j,k))*dxi
             dcdxm = (C(i,j,k)-C(im,j,k))*dxi
             dcdyp = (C(i,jp,k)-C(i,j,k))*dyi
             dcdym = (C(i,j,k)-C(i,jm,k))*dyi
             dcdzp = (C(i,j,kp)-C(i,j,k))*dzi
             dcdzm = (C(i,j,k)-C(i,j,km))*dzi
             !
             ! Scalar balance
             !
             dcdt(i,j,k) = dxi*( -ucip + ucim ) + kappa*(dcdxp-dcdxm)*dxi + &
                           dyi*( -vcjp + vcjm ) + kappa*(dcdyp-dcdym)*dyi + &
                           dzi*( -wckp + wckm ) + kappa*(dcdzp-dcdzm)*dzi
          enddo
       enddo
    enddo
   !$omp end parallel
   !
   return
  end subroutine scalad
end module mod_scal