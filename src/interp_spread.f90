module mod_interp_spread
  use mod_param
  use mod_common
  use mod_common_mpi
  use mod_kernel
  use mod_bound
  implicit none
  private
  public eulr2lagr,lagr2eulr
contains

  subroutine eulr2lagr(ibmiter)
    implicit none
    integer, intent(in) :: ibmiter
    integer :: i,j,k,l,p
    integer :: ilow,jlow,klow, ilows,jlows,klows ,ihigh,jhigh,khigh
    real :: coorx,coory,coorz
    real :: coorxs,coorys,coorzs
    real :: coorxfp,cooryfp,coorzfp
    real :: kernelx,kernely,kernelz
    real :: kernelxs,kernelys,kernelzs
    integer :: nb,nbsend,nbrecv
    integer :: nrrequests
    integer :: arrayrequests(1:3) !3=3*1 (master might have 3 slaves)
    integer :: arraystatuses(MPI_STATUS_SIZE,1:3)
    real :: boundleftnb,boundrightnb,boundfrontnb,boundbacknb
    !real, dimension(1:nl,1:npmax) :: sumu
    character*7 :: partnr,rankpr
    integer :: idp,tag
    real :: auxu,auxv,auxw
    real :: rcoeff
    real :: isperiodx,isperiody
    real :: coeffx1,coeffy1,coeffz1, &
         coeffx2,coeffy2,coeffz2, &
         coeffx3,coeffy3,coeffz3
    real :: phix(1:3), phiy(1:3), phiz(1:3), phix_half(1:3), phiy_half(1:3), phiz_half(1:3)
    !logical :: isout
    integer :: val_k_min, val_k_max

    val_k_min = 1
    val_k_max = kmax
    !
    ! third step: perform partial integration.
    !
    !$omp workshare
    forall(l=1:nl,p=1:pmax)
       ap(p)%dudtl(l) = 0.
       ap(p)%dvdtl(l) = 0.
       ap(p)%dwdtl(l) = 0.
       !    sumu(:,p)      = 0.
    end forall
    !$omp end workshare 
    !
    call updthalos_ibm(dudtf,1,1)
    call updthalos_ibm(dvdtf,1,1)
    call updthalos_ibm(dwdtf,1,1)
    call updthalos_ibm(dudtf,2,1)
    call updthalos_ibm(dvdtf,2,1)
    call updthalos_ibm(dwdtf,2,1)
    !
    !$omp parallel default(none)             &
    !$omp&shared(ap,nla,pmax)                &
    !$omp&shared(dudtf,dvdtf,dwdtf)          &
    !$omp&shared(zstart,val_k_min,val_k_max) &
    !$omp&private(p,nbrecv,l,coorxfp,cooryfp,coorzfp,ilow,ihigh,jlow,jhigh,klow,khigh) &
    !$omp&private(i,j,k,nbsend,nb) &
    !$omp&private(coorzs,coorz,kernelzs,kernelz,coorys,coory,kernelys,kernely,coorxs,coorx,kernelxs,kernelx) &
    !$omp&private(phix,phiy,phiz,phix_half,phiy_half,phiz_half) &
    !$omp&private(ilows,jlows,klows)
    !$omp do 
    do p=1,pmax
       if (ap(p)%mslv .ne. 0) then
          do l=1,nla(p)
             ! description of the coordinates on a local frame of reference with normalized space steps, 
             ! taking into account the periodicities in the axis x and y 
             coorxfp = ap(p)%xfp(l)*dxi - floor((ap(p)%xfp(l)-0.5*dx)/lx)*(itot) - (zstart(1)-1)
             cooryfp = ap(p)%yfp(l)*dyi - floor((ap(p)%yfp(l)-0.5*dy)/ly)*(jtot) - (zstart(2)-1)
             coorzfp = ap(p)%zfp(l)*dzi
             
!!! DEBUG
!             isperiodx = 0.
!             isperiody = 0. 
!             if (ap(p)%xfp(l).lt.0.+0.5*dx) isperiodx =  1.
!             if (ap(p)%xfp(l).ge.lx+0.5*dx) isperiodx = -1.
!             if (ap(p)%yfp(l).lt.0.+0.5*dy) isperiody =  1.
!             if (ap(p)%yfp(l).ge.ly+0.5*dy) isperiody = -1.
!             if (abs(coorxfp - (ap(p)%xfp(l)+isperiodx*lx-boundleftmyid )*dxi).gt.(1.e-12)) then
!                print*,'CALCULATION MISTAKE ON COORXFP ',coorxfp,(ap(p)%xfp(l)+isperiodx*lx-boundleftmyid )*dxi
!             endif
!             
!             if (abs(cooryfp - (ap(p)%yfp(l)+isperiody*ly-boundfrontmyid)*dyi ).gt.(1.e-12)) then
!                print*,'CALCULATION MISTAKE ON COORYFP ',cooryfp,  (ap(p)%yfp(l)+isperiody*ly-boundfrontmyid)*dyi
!             endif
!!! END DEBUG             

             ! Indexes for the beginning of the stencil for coorxs, coorys,coorzs
             ! (the only coordinates for which the kernel calculation is non null) 
             ilows =  nint(coorxfp) -1 ! round up
             jlows =  nint(cooryfp) -1 ! round down
             klows =  max( nint(coorzfp) -1 , val_k_min )
             ! indexes for the beginning of the stencil for coorx, coory, coorz
             ilow = floor(coorxfp)
             jlow = floor(cooryfp)
             klow = max( min( floor(coorzfp) , val_k_max-2 ) , val_k_min )
!
              
!!!THIRD VERSION :   
             phix(1) = kernel(  1.*ilows    - coorxfp)
             phix(2) = kernel(  1.*ilows +1. - coorxfp)
             phix(3) = kernel(  1.*ilows +2. - coorxfp)

             phiy(1) = kernel(  1.*jlows    - cooryfp)
             phiy(2) = kernel(  1.*jlows +1. - cooryfp)
             phiy(3) = kernel(  1.*jlows +2. - cooryfp)

             phiz(1) = kernel(  1.*klows    - coorzfp)
             phiz(2) = kernel(  1.*klows +1. - coorzfp)
             phiz(3) = kernel(  1.*klows +2. - coorzfp)

             phix_half(1) = kernel(  1.*ilow    -.5 - coorxfp)
             phix_half(2) = kernel(  1.*ilow +1. -.5 - coorxfp)
             phix_half(3) = kernel(  1.*ilow +2. -.5 - coorxfp)

             phiy_half(1) = kernel(  1.*jlow    -.5 - cooryfp)
             phiy_half(2) = kernel(  1.*jlow +1. -.5 - cooryfp)
             phiy_half(3) = kernel(  1.*jlow +2. -.5 - cooryfp)

             phiz_half(1) = kernel(  1.*klow    -.5 - coorzfp)
             phiz_half(2) = kernel(  1.*klow +1. -.5 - coorzfp)
             phiz_half(3) = kernel(  1.*klow +2. -.5 - coorzfp)

             do k=0,2
                do j=0,2
                   do i=0,2
                      ap(p)%dudtl(l) = ap(p)%dudtl(l) + dudtf(ilows+i,jlow+j,klow+k) * phix(i+1) * phiy_half(j+1) * phiz_half(k+1)
                   enddo
                enddo
             enddo

             do k=0,2
                do j=0,2
                   do i=0,2
                      ap(p)%dvdtl(l) = ap(p)%dvdtl(l) + dvdtf(ilow+i,jlows+j,klow+k) * phix_half(i+1) * phiy(j+1) * phiz_half(k+1)
                   enddo
                enddo
             enddo
             
             do k=0,2
                do j=0,2
                   do i=0,2
                      ap(p)%dwdtl(l) = ap(p)%dwdtl(l) + dwdtf(ilow+i,jlow+j,klows+k) * phix_half(i+1) * phiy_half(j+1) * phiz(k+1)
                   enddo
                enddo
             enddo
!!!END THIRD VERSION

!!! FIRST VERSION
             !isperiodx = 0.                                             
             !isperiody = 0.                                             
             !if (ap(p)%xfp(l).lt.0.+0.5*dx) isperiodx =  1.             
             !if (ap(p)%xfp(l).ge.lx+0.5*dx) isperiodx = -1.             
             !if (ap(p)%yfp(l).lt.0.+0.5*dy) isperiody =  1.             
             !if (ap(p)%yfp(l).ge.ly+0.5*dy) isperiody = -1.             
             !      isout = .false.                                      
             !coorxfp = (ap(p)%xfp(l)+isperiodx*lx-boundleftmyid )*dxi   
             !      if( nint(coorxfp).lt.1 .or. nint(coorxfp) .gt.imax ) isout = .true.
             !cooryfp = (ap(p)%yfp(l)+isperiody*ly-boundfrontmyid)*dyi 
             !      if( nint(cooryfp).lt.1 .or. nint(cooryfp) .gt.jmax ) isout = .true.
             !      if (.not.isout) then
             !ilow  = nint( coorxfp - 1.5)
             !if( ((1.*ilow     ) - coorxfp) .lt. -1.5) ilow  = ilow + 1
             !ihigh = nint( coorxfp + 1.5)
             !if( ((1.*ihigh-0.5) - coorxfp) .gt.  1.5) ihigh = ihigh - 1
             !jlow  = nint( cooryfp - 1.5)
             !if( ((1.*jlow     ) - cooryfp) .lt. -1.5) jlow  = jlow + 1
             !jhigh = nint( cooryfp + 1.5)
             !if( ((1.*jhigh-0.5) - cooryfp) .gt.  1.5) jhigh = jhigh - 1
             !klow  = nint( coorzfp - 1.5)
             !if( ((1.*klow     ) - coorzfp) .lt. -1.5) klow  = klow + 1
             !khigh = nint( coorzfp + 1.5)
             !if( ((1.*khigh-0.5) - coorzfp) .gt.  1.5) khigh = khigh - 1
             !        if (ilow .lt. -1) ilow = -1
             !        if (ihigh .gt. i1+1) ihigh = i1+1
             !        if (jlow .lt. -1) jlow = -1
             !        if (jhigh .gt. j1+1) jhigh = j1+1
             !if (klow .lt. 1) klow = 1 
             !if (khigh .gt. kmax) khigh = kmax
 
             !do k=klow,khigh
             !   coorzs   = (1.*k)-coorzfp
             !   coorz    = coorzs-0.5           !vert.    distance in grid points
             !   kernelzs = kernel(coorzs)
             !   kernelz  = kernel(coorz)
             !   coeffx1  = kernelz
             !   coeffy1  = kernelz
             !   coeffz1  = kernelzs
             !   do j=jlow,jhigh
             !      coorys   = (1.*j)-cooryfp
             !      coory    = coorys-0.5         !spanw.   distance in grid points
             !      kernelys = kernel(coorys)
             !      kernely  = kernel(coory)
             !      coeffx2  = coeffx1*kernely
             !      coeffy2  = coeffy1*kernelys
             !      coeffz2  = coeffz1*kernely
             !      do i=ilow,ihigh
             !         coorxs   = (1.*i)-coorxfp
             !         coorx    = coorxs-0.5       !streamw. distance in grid points
             !         kernelxs = kernel(coorxs)
             !         kernelx  = kernel(coorx)
             !         coeffx3  = coeffx2*kernelxs
             !         coeffy3  = coeffy2*kernelx
             !         coeffz3  = coeffz2*kernelx
             !         ap(p)%dudtl(l) = ap(p)%dudtl(l) + dudtf(i,j,k)*coeffx3
             !         ap(p)%dvdtl(l) = ap(p)%dvdtl(l) + dvdtf(i,j,k)*coeffy3
             !         ap(p)%dwdtl(l) = ap(p)%dwdtl(l) + dwdtf(i,j,k)*coeffz3
             !         !              sumu(l,p) = sumu(l,p) + kernelxs*kernely*kernelz
             !      enddo
             !   enddo
             !enddo
!!! END FIRST VERSION                                                                                                                                                   
             !!      endif
          enddo ! do l=
       endif
    enddo
    !$omp end parallel
    !
    return
  end subroutine eulr2lagr
  !
  subroutine lagr2eulr
    implicit none
    integer :: i,j,k,l,p
    integer :: ilow,jlow,klow, ilows,jlows,klows, ihigh,jhigh,khigh
    real :: coorx,coory,coorz
    real :: coorxs,coorys,coorzs
    real :: coorxfp,cooryfp,coorzfp
    real :: kernelx,kernely,kernelz
    real :: kernelxs,kernelys,kernelzs
    real :: forcex_sc,forcey_sc,forcez_sc
    integer :: nb,nbsend,nbrecv
    integer :: nrrequests
    integer :: arrayrequests(1:3) !3=3*1 (master might have 3 slaves)
    integer :: arraystatuses(MPI_STATUS_SIZE,1:3)
    real :: boundleftnb,boundrightnb,boundfrontnb,boundbacknb
    real :: dVlagrdVeuli
    !real, dimension(1:nl,1:npmax) :: sumu
    integer :: idp,tag
    real :: isperiodx,isperiody
    real :: coeffx0,coeffy0,coeffz0, &
         coeffx1,coeffy1,coeffz1, &
         coeffx2,coeffy2,coeffz2, &
         coeffx3,coeffy3,coeffz3
    !logical :: isout
    real :: phix(1:3), phiy(1:3), phiz(1:3), phix_half(1:3), phiy_half(1:3), phiz_half(1:3)
    integer :: val_k_min, val_k_max
    real :: forceytot_all
    val_k_min = 1
    val_k_max = kmax

    !
    dVlagrdVeuli = dVlagr/dVeul
    !
    !$omp workshare
    !forcex(:,:,:) = 0.
    !forcey(:,:,:) = 0.
    !forcez(:,:,:) = 0.
    dudtf(:,:,:) = 0.
    dvdtf(:,:,:) = 0.
    dwdtf(:,:,:) = 0.
    !$omp end workshare
    coeffx0 = dVlagrdVeuli*dt
    coeffy0 = dVlagrdVeuli*dt
    coeffz0 = dVlagrdVeuli*dt
    forceytot = 0.
    !
    !$omp parallel default(none)             &
    !$omp&shared(ap,nla,pmax)                &
    !$omp&shared(dudtf,dvdtf,dwdtf)          &     
    !$omp&shared(coeffx0,coeffy0,coeffz0)    &
    !$omp&shared(zstart,val_k_min,val_k_max) &
    !$omp&private(p,nbrecv,l,coorxfp,cooryfp,coorzfp,ilow,ihigh,jlow,jhigh,klow,khigh) &
    !$omp&private(i,j,k,nbsend,nb) &
    !$omp&private(coorzs,coorz,kernelzs,kernelz,coorys,coory,kernelys,kernely,coorxs,coorx,kernelxs,kernelx,forcex_sc,forcey_sc,forcez_sc) &
    !$omp&private(phix,phiy,phiz,phix_half,phiy_half,phiz_half) &
    !$omp&private(ilows,jlows,klows)
    !$omp do 
    do p=1,pmax
       if (ap(p)%mslv .ne. 0) then
          do l=1,nla(p)
             !!!THIRD VERSION
             coorxfp = ap(p)%xfp(l)*dxi - floor((ap(p)%xfp(l)-0.5*dx)/lx)*(itot) - (zstart(1)-1)
             cooryfp = ap(p)%yfp(l)*dyi - floor((ap(p)%yfp(l)-0.5*dy)/ly)*(jtot) - (zstart(2)-1)
             coorzfp = ap(p)%zfp(l)*dzi
             
             ! Indexes for the beginning of the stencil for coorxs, coorys,coorzs 
             ! (the only coordinates for which the kernel calculation is non null)
             ilows =  nint(coorxfp) -1
             jlows =  nint(cooryfp) -1
             klows =  max( nint(coorzfp) -1 , val_k_min )
             ! indexes for the beginning of the stencil for coorx, coory, coorz 
             ilow = floor(coorxfp)
             jlow = floor(cooryfp)
             klow = max( min( floor(coorzfp) , val_k_max-2 ) , val_k_min )
!
             phix(1) = kernel(  1.*ilows    - coorxfp)
             phix(2) = kernel(  1.*ilows +1. - coorxfp)
             phix(3) = kernel(  1.*ilows +2. - coorxfp)

             phiy(1) = kernel(  1.*jlows    - cooryfp)
             phiy(2) = kernel(  1.*jlows +1. - cooryfp)
             phiy(3) = kernel(  1.*jlows +2. - cooryfp)

             phiz(1) = kernel(  1.*klows    - coorzfp)
             phiz(2) = kernel(  1.*klows +1. - coorzfp)
             phiz(3) = kernel(  1.*klows +2. - coorzfp)

             phix_half(1) = kernel(  1.*ilow    -.5 - coorxfp)
             phix_half(2) = kernel(  1.*ilow +1. -.5 - coorxfp)
             phix_half(3) = kernel(  1.*ilow +2. -.5 - coorxfp)

             phiy_half(1) = kernel(  1.*jlow    -.5 - cooryfp)
             phiy_half(2) = kernel(  1.*jlow +1. -.5 - cooryfp)
             phiy_half(3) = kernel(  1.*jlow +2. -.5 - cooryfp)

             phiz_half(1) = kernel(  1.*klow    -.5 - coorzfp)
             phiz_half(2) = kernel(  1.*klow +1. -.5 - coorzfp)
             phiz_half(3) = kernel(  1.*klow +2. -.5 - coorzfp)

             do k=0,2
                do j=0,2
                   do i=0,2
                      !$omp atomic
                      dudtf(ilows+i,jlow+j,klow+k) = dudtf(ilows+i,jlow+j,klow+k) + ap(p)%fxl(l) * coeffx0 * phix(i+1) * phiy_half(j+1) * phiz_half(k+1)
                   enddo
                enddo
             enddo

             do k=0,2
                do j=0,2
                   do i=0,2
                      !$omp atomic
                      dvdtf(ilow+i,jlows+j,klow+k) = dvdtf(ilow+i,jlows+j,klow+k) + ap(p)%fyl(l) * coeffy0 * phix_half(i+1) * phiy(j+1) * phiz_half(k+1)
                      forceytot = forceytot + ap(p)%fyl(l) * coeffy0 * phix_half(i+1) * phiy(j+1) * phiz_half(k+1)
                   enddo
                enddo
             enddo

             do k=0,2
                do j=0,2
                   do i=0,2
                      !$omp atomic
                      dwdtf(ilow+i,jlow+j,klows+k) = dwdtf(ilow+i,jlow+j,klows+k) + ap(p)%fzl(l) * coeffz0 * phix_half(i+1) * phiy_half(j+1) * phiz(k+1)
                   enddo
                enddo
             enddo

             !!!FIRST VERSION
             !isperiodx = 0.
             !isperiody = 0.
             !if (ap(p)%xfp(l).lt.0.+0.5*dx) isperiodx =  1.
             !if (ap(p)%xfp(l).ge.lx+0.5*dx) isperiodx = -1.
             !if (ap(p)%yfp(l).lt.0.+0.5*dy) isperiody =  1.
             !if (ap(p)%yfp(l).ge.ly+0.5*dy) isperiody = -1.
             !!      isout = .false.
             !coorxfp = (ap(p)%xfp(l)+isperiodx*lx-boundleftmyid )*dxi
             !!      if( nint(coorxfp).lt.1 .or. nint(coorxfp) .gt.imax ) isout = .true.
             !cooryfp = (ap(p)%yfp(l)+isperiody*ly-boundfrontmyid)*dyi
             !!      if( nint(cooryfp).lt.1 .or. nint(cooryfp) .gt.jmax ) isout = .true.
             !!      if (.not.isout) then
             !coorzfp = ap(p)%zfp(l)*dzi
             !ilow  = nint( coorxfp - 1.5)
             !if( ((1.*ilow     ) - coorxfp) .lt. -1.5) ilow  = ilow + 1
             !ihigh = nint( coorxfp + 1.5)
             !if( ((1.*ihigh-0.5) - coorxfp) .gt.  1.5) ihigh = ihigh - 1
             !jlow  = nint( cooryfp - 1.5)
             !if( ((1.*jlow     ) - cooryfp) .lt. -1.5) jlow  = jlow + 1
             !jhigh = nint( cooryfp + 1.5)
             !if( ((1.*jhigh-0.5) - cooryfp) .gt.  1.5) jhigh = jhigh - 1
             !klow  = nint( coorzfp - 1.5)
             !if( ((1.*klow     ) - coorzfp) .lt. -1.5) klow  = klow + 1
             !khigh = nint( coorzfp + 1.5)
             !if( ((1.*khigh-0.5) - coorzfp) .gt.  1.5) khigh = khigh - 1
             !!        if (ilow .lt. -1) ilow = -1
             !!        if (ihigh .gt. i1+1) ihigh = i1+1
             !!        if (jlow .lt. -1) jlow = -1 
             !!        if (jhigh .gt. j1+1) jhigh = j1+1 
             !if (klow .lt. 1) klow = 1
             !if (khigh .gt. kmax) khigh = kmax
            ! 
            ! do k=klow,khigh
            !    coorzs   = (1.*k)-coorzfp
            !    coorz    = coorzs-0.5           !vert.    distance in grid points
            !    kernelzs = kernel(coorzs)
            !    kernelz  = kernel(coorz)
            !    coeffx1   = coeffx0*kernelz
            !    coeffy1   = coeffy0*kernelz
            !    coeffz1   = coeffz0*kernelzs
            !    do j=jlow,jhigh
            !       coorys   = (1.*j)-cooryfp
            !       coory    = coorys-0.5         !spanw.   distance in grid points
            !       kernelys = kernel(coorys)
            !       kernely  = kernel(coory)
            !       coeffx2   = coeffx1*kernely
            !       coeffy2   = coeffy1*kernelys
            !       coeffz2   = coeffz1*kernely
            !       do i=ilow,ihigh
            !          coorxs   = (1.*i)-coorxfp
            !          coorx    = coorxs-0.5       !streamw. distance in grid points
            !          kernelxs = kernel(coorxs)
            !          kernelx  = kernel(coorx)
            !          coeffx3  = coeffx2*kernelxs
            !          coeffy3  = coeffy2*kernelx
            !          coeffz3  = coeffz2*kernelx
            !          !  !$omp atomic
            !          !              forcex(i,j,k) = forcex(i,j,k) + ap(p)%fxl(l)*kernelxs*kernely*kernelz*dVlagrdVeuli
            !          !  !$omp atomic
            !          !              forcey(i,j,k) = forcey(i,j,k) + ap(p)%fyl(l)*kernelx*kernelys*kernelz*dVlagrdVeuli
            !          !  !$omp atomic
            !          !              forcez(i,j,k) = forcez(i,j,k) + ap(p)%fzl(l)*kernelx*kernely*kernelzs*dVlagrdVeuli
            !          !$omp atomic
            !             !              forcex(i,j,k) = forcex(i,j,k) + ap(p)%fxl(l)*kernelxs*kernely*kernelz*dVlagrdVeuli
            !          !  !$omp atomic
            !          !              forcey(i,j,k) = forcey(i,j,k) + ap(p)%fyl(l)*kernelx*kernelys*kernelz*dVlagrdVeuli
            !          !  !$omp atomic
            !          !              forcez(i,j,k) = forcez(i,j,k) + ap(p)%fzl(l)*kernelx*kernely*kernelzs*dVlagrdVeuli
            !          !$omp atomic
            !          dudtf(i,j,k) = dudtf(i,j,k) + ap(p)%fxl(l)*coeffx3
            !          !$omp atomic
            !          dvdtf(i,j,k) = dvdtf(i,j,k) + ap(p)%fyl(l)*coeffy3
            !          !$omp atomic
            !          dwdtf(i,j,k) = dwdtf(i,j,k) + ap(p)%fzl(l)*coeffz3
            !       enddo
            !    enddo
            ! enddo
            ! !       endif
          enddo !do l=
       endif
    enddo
    !$omp end parallel
    
    !
    call updthalos_ibm(dudtf ,1,2)
    call updthalos_ibm(dvdtf ,1,2)
    call updthalos_ibm(dwdtf ,1,2)
    call updthalos_ibm(dudtf ,2,2)
    call updthalos_ibm(dvdtf ,2,2)
    call updthalos_ibm(dwdtf ,2,2)
    !
    call mpi_allreduce(forceytot,forceytot_all,1,mpi_real8,mpi_sum,comm_cart,error)
    forceytot = forceytot_all/dt/(1.*itot*jtot*ktot) + wallshearnew
    if(iniu.ne.'COU') then
      do k=1,kmax
        do j=1,jmax
          do i=1,imax
            dudtf(i,j,k) = dudtf(i,j,k) + dudtold(i,j,k)
            dvdtf(i,j,k) = dvdtf(i,j,k) + dvdtold(i,j,k) ! + (-forceytot)*dt ! free falling, no need to balance lift force in y-direction (SS)
            dwdtf(i,j,k) = dwdtf(i,j,k) + dwdtold(i,j,k)
          enddo
        enddo
      enddo
    endif
    !call updthalos_ibm(forcex,1,2)
    !call updthalos_ibm(forcey,1,2)
    !call updthalos_ibm(forcez,1,2)
    !call updthalos_ibm(forcex,2,2)
    !call updthalos_ibm(forcey,2,2)
    !call updthalos_ibm(forcez,2,2)
    !
!    !$omp workshare
!    dvdtf(1:imax,1:jmax,1:kmax) = dvdtf(1:imax,1:jmax,1:kmax) + 1.0*(bulk_v_sup-v_bulk)
!    !$omp end workshare
    !
    return
  end subroutine lagr2eulr
  !
end module mod_interp_spread
