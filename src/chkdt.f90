module mod_chkdt
  implicit none
  private
  public chkdt
contains
  subroutine chkdt(dtmax)
    use mpi
    use mod_param
    use mod_common
    use mod_common_mpi
    implicit none
    !
    ! calculates the timestep dt based on stability conditions
    ! for RK3 scheme (see: P. Wesseling, 'Principles of computational
    ! fluid dynamics').
    !
    real, intent(inout) :: dtmax
    real :: velo,dxi2
    real :: u2,v2,w2
    real :: dt1,dt2
    real :: dtmaxold,dtmax_all(2,1)
    real :: dtmaxu(2,1)
    real :: dtmaxv(2,1)
    real :: dtmaxw(2,1)
    real  :: flagmaxu,flagmaxv,flagmaxw
    real  :: dummy
    integer :: counter1,counter2
    integer :: counter1_all,counter2_all
    integer :: counter1u,counter2u
    integer :: counter1v,counter2v
    integer :: counter1w,counter2w
    integer :: maxi,maxj,maxk
    integer :: maxiu,maxju,maxku
    integer :: maxiv,maxjv,maxkv
    integer :: maxiw,maxjw,maxkw
    integer :: flag,flagmax
    integer :: i,j,k
    !
    ! streamwise velocity
    !
    dtmaxu(1,1) = 899999.
    dtmaxu(2,1) = myid*1.
    dtmaxold    = dtmaxu(1,1)
    counter1    = 0
    counter2    = 0
    maxi        = 1000
    maxj        = 1000
    maxk        = 1000
    flag        = 0
    flagmax     = 0
        
    !!$omp parallel default(none) &
    !!$omp&shared(unew,vnew,wnew) &
    !!$omp&private(i,j,k,dxi2,dt1,u2,v2,w2,velo,dt2,flag,counter1,counter2,dummy) &
    !!$omp&reducion(min:)
    !!$omp do
    do k=1,kmax
       do j=1,jmax
          do i=1,imax
             ! diffusion
             dxi2  = dxi**2 + dyi**2 + dzi**2
             dt1   = 1.65/( 4.*visc*dxi2 )  ! 2.5127/( 4.*visc*dxi2 )
             ! advection
             u2 = ( unew(i,j,k) )**2.
             v2 = ( 0.25*(vnew(i,j,k)+vnew(i,j-1,k)+vnew(i+1,j,k)+vnew(i+1,j-1,k)) )**2.
             w2 = ( 0.25*(wnew(i,j,k)+wnew(i,j,k-1)+wnew(i+1,j,k)+wnew(i+1,j,k-1)) )**2.
             velo = sqrt(u2+v2+w2)
             if (abs(1.*velo) .gt. 1.e-12) then
                dt2   = ( sqrt(1.5) )/( velo*sqrt(dxi2) ) ! Isotropic condition (SS)
                dt2   = sqrt(3.)/( sqrt(u2)*dxi + sqrt(v2)*dyi + sqrt(w2)*dzi ) ! Anisotropy correction (SS)
             else
                dt2 = 99999.
             endif
             if (dt1 .lt. dt2) then
                flag = 1
                counter1 = counter1 + 1
             else
                flag = 2
                counter2 = counter2 + 1
             endif
             dummy       = dtmaxu(1,1)
             dtmaxu(1,1) = min(dummy,dt1,dt2)
             if (dtmaxu(1,1) .lt. dtmaxold) then
                maxi     = i
                maxj     = j
                maxk     = k
                flagmax  = flag
                dtmaxold = dtmaxu(1,1)
             endif
          enddo
       enddo
    enddo

    !!$omp end parallel
    call mpi_allreduce(dtmaxu,dtmax_all,1,mpi_2double_precision,mpi_minloc,comm_cart,error)
    call mpi_allreduce(counter1,counter1_all,1,mpi_integer,mpi_sum,comm_cart,error)
    call mpi_allreduce(counter2,counter2_all,1,mpi_integer,mpi_sum,comm_cart,error)
    dtmaxu(1,1) = dtmax_all(1,1)
    dtmaxu(2,1) = dtmax_all(2,1)
    maxiu       = maxi
    maxju       = maxj
    maxku       = maxk
    counter1u   = counter1_all
    counter2u   = counter2_all
    flagmaxu    = flagmax
    !
    ! spanwise velocity
    !
    dtmaxv(1,1) = 899999.
    dtmaxv(2,1) = myid*1.
    dtmaxold    = dtmaxv(1,1)
    counter1    = 0
    counter2    = 0
    maxi        = 1000
    maxj        = 1000
    maxk        = 1000
    flag        = 0
    flagmax     = 0
    do k=1,kmax
       do j=1,jmax
          do i=1,imax
             ! diffusion
             dxi2  = dxi**2 + dyi**2 + dzi**2
             dt1   = 1.65/( 4.*visc*dxi2 ) ! 2.5127/( 4.*visc*dxi2 )
             ! advection
             u2 = ( 0.25*(unew(i,j,k)+unew(i,j+1,k)+unew(i-1,j+1,k)+unew(i-1,j,k)) )**2.
             v2 = ( vnew(i,j,k) )**2.
             w2 = ( 0.25*(wnew(i,j,k)+wnew(i,j+1,k)+wnew(i,j+1,k-1)+wnew(i,j,k-1)) )**2.
             velo = sqrt(u2+v2+w2)
             if (abs(1.*velo) .gt. 1.e-12) then
                dt2   = ( sqrt(1.5) )/( velo*sqrt(dxi2) )
                dt2   = sqrt(3.)/( sqrt(u2)*dxi + sqrt(v2)*dyi + sqrt(w2)*dzi )
             else
                dt2 = 99999.
             endif
             if (dt1 .lt. dt2) then
                flag = 1
                counter1 = counter1 + 1
             else
                flag = 2
                counter2 = counter2 + 1
             endif
             dummy       = dtmaxv(1,1)
             dtmaxv(1,1) = min(dummy,dt1,dt2)
             if (dtmaxv(1,1) .lt. dtmaxold) then
                maxi     = i
                maxj     = j
                maxk     = k
                flagmax  = flag
                dtmaxold = dtmaxv(1,1)
             endif
          enddo
       enddo
    enddo
    call mpi_allreduce(dtmaxv,dtmax_all,1,mpi_2double_precision,mpi_minloc,comm_cart,error)
    call mpi_allreduce(counter1,counter1_all,1,mpi_integer,mpi_sum,comm_cart,error)
    call mpi_allreduce(counter2,counter2_all,1,mpi_integer,mpi_sum,comm_cart,error)
    dtmaxv(1,1) = dtmax_all(1,1)
    dtmaxv(2,1) = dtmax_all(2,1)
    maxiv       = maxi
    maxjv       = maxj
    maxkv       = maxk
    counter1v   = counter1_all
    counter2v   = counter2_all
    flagmaxv    = flagmax
    !
    ! wall-normal velocity
    !
    dtmaxw(1,1) = 899999.
    dtmaxw(2,1) = myid*1.
    dtmaxold    = dtmaxw(1,1)
    counter1    = 0
    counter2    = 0
    maxi        = 1000
    maxj        = 1000
    maxk        = 1000
    flag        = 0
    flagmax     = 0
    do k=1,kmax
       do j=1,jmax
          do i=1,imax
             ! diffusion
             dxi2  = dxi**2 + dyi**2 + dzi**2
             dt1   = 1.65/( 4.*visc*dxi2 ) ! 2.5127/( 4.*visc*dxi2 )
             ! advection
             u2 = ( 0.25*(unew(i,j,k)+unew(i-1,j,k)+unew(i-1,j,k+1)+unew(i,j,k+1)) )**2.
             v2 = ( 0.25*(vnew(i,j,k)+vnew(i,j-1,k)+vnew(i,j-1,k+1)+vnew(i,j,k+1)) )**2.
             w2 = ( wnew(i,j,k) )**2.
             velo = sqrt(u2+v2+w2)
             if (abs(1.*velo) .gt. 1.e-12) then
                dt2   = ( sqrt(1.5) )/( velo*sqrt(dxi2) )
                dt2   = sqrt(3.)/( sqrt(u2)*dxi + sqrt(v2)*dyi + sqrt(w2)*dzi )
             else
                dt2 = 99999.
             endif
             if (dt1 .lt. dt2) then
                flag = 1
                counter1 = counter1 + 1
             else
                flag = 2
                counter2 = counter2 + 1
             endif
             dummy       = dtmaxw(1,1)
             dtmaxw(1,1) = min(dummy,dt1,dt2)
             if (dtmaxw(1,1) .lt. dtmaxold) then
                maxi     = i
                maxj     = j
                maxk     = k
                flagmax  = flag
                dtmaxold = dtmaxw(1,1)
             endif
          enddo
       enddo
    enddo
    call mpi_allreduce(dtmaxw,dtmax_all,1,mpi_2double_precision,mpi_minloc,comm_cart,error)
    call mpi_allreduce(counter1,counter1_all,1,mpi_integer,mpi_sum,comm_cart,error)
    call mpi_allreduce(counter2,counter2_all,1,mpi_integer,mpi_sum,comm_cart,error)
    dtmaxw(1,1) = dtmax_all(1,1)
    dtmaxw(2,1) = dtmax_all(2,1)
    maxiw       = maxi
    maxjw       = maxj
    maxkw       = maxk
    counter1w   = counter1_all
    counter2w   = counter2_all
    flagmaxw    = flagmax
    !
    ! determine maximum time step
    !
    dtmax = min(dtmaxu(1,1), dtmaxv(1,1), dtmaxw(1,1))
    !
    ! x-momentum
    !
    if ( dtmaxu(1,1) .eq. dtmax ) then
       if (myid .eq. nint(dtmaxu(2,1))) then
          write(6,*) 'STABILITY FOR STREAMWISE MOMENTUM EQUATION'
          write(6,*) 'Maximum allowed delta t, i, j, k = ', dtmaxu(1,1), maxiu+coords(1)*itot/dims(1), &
               maxju+coords(2)*jtot/dims(2), maxku 
          write(6,*) '%1, %2, flag = ', &
               100*real(counter1u)/(itot*jtot*kmax), &
               100*real(counter2u)/(itot*jtot*kmax),flagmaxu
       endif
    endif
    !
    ! y-momentum
    !
    if ( dtmaxv(1,1) .eq. dtmax ) then
       if (myid .eq. nint(dtmaxv(2,1))) then
          write(6,*) 'STABILITY FOR SPANWISE MOMENTUM EQUATION'
          write(6,*) 'Maximum allowed delta t, i, j, k = ', dtmaxv(1,1), maxiv+coords(1)*itot/dims(1), &
               maxjv+coords(2)*jtot/dims(2), maxkv
          write(6,*) '%1, %2, flag = ', &
               100*real(counter1v)/(itot*jtot*kmax), &
               100*real(counter2v)/(itot*jtot*kmax),flagmaxv
       endif
    endif
    !
    ! z-momentum
    !
    if ( dtmaxw(1,1) .eq. dtmax ) then
       if (myid .eq. nint(dtmaxw(2,1))) then
          write(6,*) 'STABILITY FOR WALL-NORMAL MOMENTUM EQUATION'
          write(6,*) 'Maximum allowed delta t, i, j, k = ', dtmaxw(1,1), maxiw+coords(1)*itot/dims(1), &
               maxjw+coords(2)*jtot/dims(2), maxkw
          write(6,*) '%1, %2, flag = ', &
               100*real(counter1w)/(itot*jtot*kmax), &
               100*real(counter2w)/(itot*jtot*kmax),flagmaxw
       endif
    endif
    !
    return
  end subroutine chkdt
  !
end module mod_chkdt
