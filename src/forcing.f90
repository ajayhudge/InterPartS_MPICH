module mod_forcing
  use mod_common
  use mod_common_mpi
  implicit none
  private
  public complagrforces,updtlagrforces,updtintermediatevel
contains
  !
  subroutine complagrforces
    implicit none
    integer :: l,p
    real :: dti
    type pneighbor
       real :: intu,intv,intw, &
               intomx,intomy,intomz
    end type pneighbor
    type(pneighbor), dimension(0:8,1:npmax) :: anb
    integer :: nb,nbsend,nbrecv
    integer :: nrrequests
    integer :: arrayrequests(1:3) ! 3=3*1 (master might have 3 slaves)
    integer :: arraystatuses(MPI_STATUS_SIZE,1:3)
    integer :: idp,tag
    real :: auxu,auxv,auxw,auxomx,auxomy,auxomz
    real :: isperiodx,isperiody
    real :: coorxfp,cooryfp,coorzfp
    logical :: isout

    dti = 1./dt
    !
    !$omp parallel default(none)         &
    !$omp&shared(ap,pmax,nla,dvlagr,dti) &
    !$omp&private(p,l)
    !$omp do
    do p=1,pmax
       if(ap(p)%mslv.ne.0) then
          !ap(p)%fxltot_old  = ap(p)%fxltot ***ask wp why this was here before***
          !ap(p)%fyltot_old  = ap(p)%fyltot
          !ap(p)%fzltot_old  = ap(p)%fzltot
          !ap(p)%torqxltotold  = ap(p)%torqxltot
          !ap(p)%torqyltotold  = ap(p)%torqyltot
          !ap(p)%torqzltotold  = ap(p)%torqzltot
          !ap(p)%torqthetaold = ap(p)%torqtheta
          ap(p)%fxltot  = 0.
          ap(p)%fyltot  = 0.
          ap(p)%fzltot  = 0.
          ap(p)%torqxltot  = 0.
          ap(p)%torqyltot  = 0.
          ap(p)%torqzltot  = 0.
          do l=1,NLa(p)
!             isperiodx = 0.
!             isperiody = 0.
!             if (ap(p)%xfp(l).lt.0.+0.5*dx) isperiodx =  1.
!             if (ap(p)%xfp(l).ge.lx+0.5*dx) isperiodx = -1.
!             if (ap(p)%yfp(l).lt.0.+0.5*dy) isperiody =  1.
!             if (ap(p)%yfp(l).ge.ly+0.5*dy) isperiody = -1.
!             isout = .false.
!             coorxfp = (ap(p)%xfp(l)+isperiodx*lx-boundleftmyid )*dxi
!             if( nint(coorxfp).lt.1 .or. nint(coorxfp) .gt.imax ) isout = .true.
!             cooryfp = (ap(p)%yfp(l)+isperiody*ly-boundfrontmyid)*dyi
!             if( nint(cooryfp).lt.1 .or. nint(cooryfp) .gt.jmax ) isout = .true.
!             if (.not.isout) then
!                coorzfp =  ap(p)%zfp(l)*dzi
                ap(p)%fxl(l) = (ap(p)%ul(l)-ap(p)%dudtl(l))*dti
                ap(p)%fyl(l) = (ap(p)%vl(l)-ap(p)%dvdtl(l))*dti
                ap(p)%fzl(l) = (ap(p)%wl(l)-ap(p)%dwdtl(l))*dti
                ap(p)%fxltot = ap(p)%fxltot + ap(p)%fxl(l)
                ap(p)%fyltot = ap(p)%fyltot + ap(p)%fyl(l)
                ap(p)%fzltot = ap(p)%fzltot + ap(p)%fzl(l)
                ap(p)%torqxltot = ap(p)%torqxltot + (ap(p)%yfp(l)-ap(p)%y)*ap(p)%fzl(l) - &
                                                    (ap(p)%zfp(l)-ap(p)%z)*ap(p)%fyl(l)
                ap(p)%torqyltot = ap(p)%torqyltot + (ap(p)%zfp(l)-ap(p)%z)*ap(p)%fxl(l) - &
                                                    (ap(p)%xfp(l)-ap(p)%x)*ap(p)%fzl(l)
                ap(p)%torqzltot = ap(p)%torqzltot + (ap(p)%xfp(l)-ap(p)%x)*ap(p)%fyl(l) - &
                                                    (ap(p)%yfp(l)-ap(p)%y)*ap(p)%fxl(l)
!             endif
          enddo
          ap(p)%fxltot = ap(p)%fxltot*dVlagr
          ap(p)%fyltot = ap(p)%fyltot*dVlagr
          ap(p)%fzltot = ap(p)%fzltot*dVlagr
          ap(p)%torqxltot = ap(p)%torqxltot*dVlagr
          ap(p)%torqyltot = ap(p)%torqyltot*dVlagr
          ap(p)%torqzltot = ap(p)%torqzltot*dVlagr
          ! Torque working on angle theta:
          ! phi = 0      --> torqtheta =  torqyltot
          ! phi = pi/2   --> torqtheta = -torqxltot
          ! phi = pi     --> torqtheta = -torqyltot
          ! phi = 3*pi/2 --> torqtheta =  torqxltot
          ap(p)%torqtheta = (ap(p)%torqyltot*cos(ap(p)%phi)) - &
               (ap(p)%torqxltot*sin(ap(p)%phi))
          ! phi is value of phic at time step n (weak coupling)
       endif
    enddo
    !$omp end parallel
    !
    !$omp workshare
    forall(p=1:npmax,nb=1:8)
       anb(nb,p)%intu   = 0.
       anb(nb,p)%intv   = 0.
       anb(nb,p)%intw   = 0.
       anb(nb,p)%intomx = 0.
       anb(nb,p)%intomy = 0.
       anb(nb,p)%intomz = 0.
    end forall
    !$omp end workshare
    !!$omp workshare
    anb(0,1:npmax)%intu   = ap(1:npmax)%fxltot 
    anb(0,1:npmax)%intv   = ap(1:npmax)%fyltot
    anb(0,1:npmax)%intw   = ap(1:npmax)%fzltot
    anb(0,1:npmax)%intomx = ap(1:npmax)%torqxltot
    anb(0,1:npmax)%intomy = ap(1:npmax)%torqyltot
    anb(0,1:npmax)%intomz = ap(1:npmax)%torqzltot
    !!$omp end workshare
    !
    do p=1,pmax
       nrrequests = 0
       do nb=1,8
          nbsend = nb    ! rank of process which sends data ('data is received from neighbor nbsend')
          idp = abs(ap(p)%mslv)
          !tag = idp*10+nbsend
          tag = idp*10+nbsend-idp*10
          nbrecv  = nb+4  ! rank of process which receives data ('data is send to neighbor nbrecv')
          if (nbrecv .gt. 8) nbrecv = nbrecv - 8
          if (ap(p)%mslv .gt. 0) then
             ! myid is master of particle ap(p)%mslv
             if (ap(p)%nb(nbsend) .eq. 1) then
                ! neighbor(nbsend) is rank of slave for particle ap(p)%mslv
                if ( neighbor(nbsend) .ne. myid ) then
                   nrrequests = nrrequests + 1
                   call MPI_IRECV(anb(nbrecv,p)%intu,6,MPI_REAL8,neighbor(nbsend), &
                        tag,comm_cart,arrayrequests((nrrequests-1) + 1),error)
                   ! recv intu,intv,intw (...) -> 6 contiguous info
                   ! (see definition of type pneighbor in the begining of the subroutine)
                endif
             endif
          endif
          if (ap(p)%mslv .lt. 0) then
             ! myid is slave of particle -ap(p)%mslv
             if (ap(p)%nb(nbrecv) .eq. 1) then
                ! neighbor(nbrecv) is rank of master of particle -ap(p)%mslv
                nrrequests = nrrequests + 1
                call MPI_ISEND(anb(0     ,p)%intu,6,MPI_REAL8,neighbor(nbrecv), &
                     tag,comm_cart,arrayrequests((nrrequests-1) + 1),error)
                ! send intu,intv,intw (...) -> 6 contiguous info
                ! (see definition of type pneighbor in the begining of the subroutine)
             endif
          endif
       enddo ! do nb=
       call MPI_WAITALL(nrrequests,arrayrequests,arraystatuses,error)
    enddo
    !$omp parallel default(none) &
    !$omp&shared(ap,anb,pmax)    &
    !$omp&private(p,nb) reduction(+:auxu,auxv,auxw,auxomx,auxomy,auxomz)
    !$omp do 
    do p=1,pmax
       ap(p)%fxltot    = 0.
       ap(p)%fyltot    = 0.
       ap(p)%fzltot    = 0.
       ap(p)%torqxltot = 0.
       ap(p)%torqyltot = 0.
       ap(p)%torqzltot = 0.
       if (ap(p)%mslv .gt. 0) then
          auxu   = 0.
          auxv   = 0.
          auxw   = 0.
          auxomx = 0.
          auxomy = 0.
          auxomz = 0.
          do nb=0,8
             auxu   = auxu   + anb(nb,p)%intu
             auxv   = auxv   + anb(nb,p)%intv
             auxw   = auxw   + anb(nb,p)%intw
             auxomx = auxomx + anb(nb,p)%intomx
             auxomy = auxomy + anb(nb,p)%intomy
             auxomz = auxomz + anb(nb,p)%intomz
          enddo
          ap(p)%fxltot    = auxu
          ap(p)%fyltot    = auxv
          ap(p)%fzltot    = auxw
          ap(p)%torqxltot = auxomx
          ap(p)%torqyltot = auxomy
          ap(p)%torqzltot = auxomz
       endif
    enddo
    !$omp end parallel
    !
    return
  end subroutine complagrforces
  !
  subroutine updtlagrforces(maxerror_all,maxavererror_all)
    implicit none
    integer l,p
    real :: maxavererror,maxerror
    real,intent(out) :: maxavererror_all,maxerror_all
    real :: relerror,avererror,lagvel,averlagvel,averlagvel_all
    real :: errdistr
    real :: dti
    type pneighbor
       real :: intu,intv,intw, &
            intomx,intomy,intomz
    end type pneighbor
    type(pneighbor), dimension(0:8,1:npmax) :: anb
    integer :: i,nb,nbsend,nbrecv
    integer :: nrrequests
    integer :: arrayrequests(1:3) ! 3=3*1 (master might have 3 slaves)
    integer :: arraystatuses(MPI_STATUS_SIZE,1:3)
    integer :: idp,tag
    real :: auxu,auxv,auxw,auxomx,auxomy,auxomz
    real :: isperiodx,isperiody
    real :: coorxfp,cooryfp,coorzfp
    logical :: isout
    !
    dti = 1./dt
    !
    maxerror = 0.
    maxavererror = 0.
    !$omp parallel default(none)         &
    !$omp&shared(ap,pmax,nla,dvlagr,dti) &
    !$omp& private(p,l,lagvel,relerror,avererror,averlagvel) &
    !$omp& reduction(max:maxerror,maxavererror)
    !$omp do
    do p=1,pmax
       if(ap(p)%mslv.ne.0) then
          ap(p)%fxltot = 0.
          ap(p)%fyltot = 0.
          ap(p)%fzltot = 0.
          ap(p)%torqxltot = 0.
          ap(p)%torqyltot = 0.
          ap(p)%torqzltot = 0.
          avererror  = 0.
          averlagvel = 0.
          do l=1,NLa(p)
!             isperiodx = 0.
!             isperiody = 0.
!             if (ap(p)%xfp(l).lt.0.+0.5*dx) isperiodx =  1.
!             if (ap(p)%xfp(l).ge.lx+0.5*dx) isperiodx = -1.
!             if (ap(p)%yfp(l).lt.0.+0.5*dy) isperiody =  1.
!             if (ap(p)%yfp(l).ge.ly+0.5*dy) isperiody = -1.
!             isout = .false.
!             coorxfp = (ap(p)%xfp(l)+isperiodx*lx-boundleftmyid )*dxi
!             if( nint(coorxfp).lt.1 .or. nint(coorxfp) .gt.imax ) isout = .true.
!             cooryfp = (ap(p)%yfp(l)+isperiody*ly-boundfrontmyid)*dyi
!             if( nint(cooryfp).lt.1 .or. nint(cooryfp) .gt.jmax ) isout = .true.
!             if (.not.isout) then
!                coorzfp =  ap(p)%zfp(l)*dzi
                !      if (abs(ap(p)%mslv) .eq. 1) then
                !        lagvel        = sqrt(ap(p)%ul(l)**2 + ap(p)%vl(l)**2 + ap(p)%wl(l)**2)
                !        errdistr      = sqrt(ap(p)%dudtl(l)**2 + ap(p)%dvdtl(l)**2 + ap(p)%dwdtl(l)**2) - &
                !                        lagvel
                !        avererror     = avererror + errdistr
                !        averlagvel    = averlagvel + lagvel
                !        relerror      = 100.*errdistr/(1.e-12 + lagvel)
                !        !if ( maxerror .lt. abs(relerror)) maxerror = abs(relerror)
                !        maxerror=max(maxerror,abs(relerror))
                !      endif
                ap(p)%fxl(l) = ap(p)%fxl(l) + (ap(p)%ul(l)-ap(p)%dudtl(l))*dti ! update force
                ap(p)%fyl(l) = ap(p)%fyl(l) + (ap(p)%vl(l)-ap(p)%dvdtl(l))*dti ! update force
                ap(p)%fzl(l) = ap(p)%fzl(l) + (ap(p)%wl(l)-ap(p)%dwdtl(l))*dti ! update force
                ap(p)%fxltot = ap(p)%fxltot + ap(p)%fxl(l)
                ap(p)%fyltot = ap(p)%fyltot + ap(p)%fyl(l)
                ap(p)%fzltot = ap(p)%fzltot + ap(p)%fzl(l)
                ap(p)%torqxltot = ap(p)%torqxltot + (ap(p)%yfp(l)-ap(p)%y)*ap(p)%fzl(l) - &
                                                    (ap(p)%zfp(l)-ap(p)%z)*ap(p)%fyl(l)
                ap(p)%torqyltot = ap(p)%torqyltot + (ap(p)%zfp(l)-ap(p)%z)*ap(p)%fxl(l) - &
                                                    (ap(p)%xfp(l)-ap(p)%x)*ap(p)%fzl(l)
                ap(p)%torqzltot = ap(p)%torqzltot + (ap(p)%xfp(l)-ap(p)%x)*ap(p)%fyl(l) - &
                                                    (ap(p)%yfp(l)-ap(p)%y)*ap(p)%fxl(l)
!             endif
          enddo
          !  if (abs(ap(p)%mslv) .eq. 1) then
          !    averlagvel = averlagvel/(1.*NL)
          !    avererror  = 100.*(avererror/(1.*NL))/(1.e-12+averlagvel) ! incorrect now!
          !  endif
          !!
          !  maxavererror=max(maxavererror,abs(avererror))
          ap(p)%fxltot = ap(p)%fxltot*dVlagr
          ap(p)%fyltot = ap(p)%fyltot*dVlagr
          ap(p)%fzltot = ap(p)%fzltot*dVlagr
          ap(p)%torqxltot = ap(p)%torqxltot*dVlagr
          ap(p)%torqyltot = ap(p)%torqyltot*dVlagr
          ap(p)%torqzltot = ap(p)%torqzltot*dVlagr
          ! Torque working on angle theta:
          !   phi = 0      --> torqtheta =  torqyltot
          !   phi = pi/2   --> torqtheta = -torqxltot
          !   phi = pi     --> torqtheta = -torqyltot
          !   phi = 3*pi/2 --> torqtheta =  torqxltot
          ap(p)%torqtheta = (ap(p)%torqyltot*cos(ap(p)%phi)) - &
                            (ap(p)%torqxltot*sin(ap(p)%phi))
          ! phi is value of phic at time step n for which torqyltot and torqxltot are computed!
       endif
    enddo
    !$omp end parallel 
    !
    !call mpi_allreduce(maxavererror,maxavererror_all,1,mpi_real8,mpi_max,comm_cart,error) ! incorrect now!
    !call mpi_allreduce(maxerror,maxerror_all,1,mpi_real8,mpi_max,comm_cart,error)
    maxavererror_all = 0.
    maxerror_all     = 0.
    !
    !$omp workshare
    forall(p=1:npmax,nb=1:8)
       anb(nb,p)%intu   = 0.
       anb(nb,p)%intv   = 0.
       anb(nb,p)%intw   = 0.
       anb(nb,p)%intomx = 0.
       anb(nb,p)%intomy = 0.
       anb(nb,p)%intomz = 0.
    end forall
    !$omp  end workshare
    !!$omp workshare
    anb(0,1:npmax)%intu   = ap(1:npmax)%fxltot 
    anb(0,1:npmax)%intv   = ap(1:npmax)%fyltot
    anb(0,1:npmax)%intw   = ap(1:npmax)%fzltot
    anb(0,1:npmax)%intomx = ap(1:npmax)%torqxltot
    anb(0,1:npmax)%intomy = ap(1:npmax)%torqyltot
    anb(0,1:npmax)%intomz = ap(1:npmax)%torqzltot
    !!$omp  end workshare
    !
    do p=1,pmax
       nrrequests = 0
       do nb=1,8
          nbsend = nb    ! rank of process which sends data ('data is received from neighbor nbsend')
          idp = abs(ap(p)%mslv)
          !tag = idp*10+nbsend
          tag = idp*10+nbsend-idp*10
          nbrecv  = nb+4  ! rank of process which receives data ('data is send to neighbor nbrecv')
          if (nbrecv .gt. 8) nbrecv = nbrecv - 8
          if (ap(p)%mslv .gt. 0) then
             ! myid is master of particle ap(p)%mslv
             if (ap(p)%nb(nbsend) .eq. 1) then
                ! neighbor(nbsend) is rank of slave for particle ap(p)%mslv
                if ( neighbor(nbsend) .ne. myid ) then
                   nrrequests = nrrequests + 1
                   call MPI_IRECV(anb(nbrecv,p)%intu,6,MPI_REAL8,neighbor(nbsend), &
                        tag,comm_cart,arrayrequests((nrrequests-1) + 1),error)
                   ! recv intu,intv,intw (...) -> 6 contiguous info
                   ! (see definition of type pneighbor in the begining of the subroutine)
                endif
             endif
          endif
          if (ap(p)%mslv .lt. 0) then
             ! myid is slave of particle -ap(p)%mslv
             if (ap(p)%nb(nbrecv) .eq. 1) then
                ! neighbor(nbrecv) is rank of master of particle -ap(p)%mslv
                nrrequests = nrrequests + 1
                call MPI_ISEND(anb(0     ,p)%intu,6,MPI_REAL8,neighbor(nbrecv), &
                     tag,comm_cart,arrayrequests((nrrequests-1) + 1),error)
                ! send intu,intv,intw (...) -> 6 contiguous info
                ! (see definition of type pneighbor in the begining of the subroutine)
             endif
          endif
       enddo ! do nb=
       call MPI_WAITALL(nrrequests,arrayrequests,arraystatuses,error)
    enddo
    !
    !$omp parallel default(none) &
    !$omp&shared(ap,anb,pmax)    &
    !$omp&private(p,nb) reduction(+:auxu,auxv,auxw,auxomx,auxomy,auxomz)
    !$omp do 
    do p=1,pmax
       ap(p)%fxltot    = 0.
       ap(p)%fyltot    = 0.
       ap(p)%fzltot    = 0.
       ap(p)%torqxltot = 0.
       ap(p)%torqyltot = 0.
       ap(p)%torqzltot = 0.
       if (ap(p)%mslv .gt. 0) then
          auxu   = 0.
          auxv   = 0.
          auxw   = 0.
          auxomx = 0.
          auxomy = 0.
          auxomz = 0.
          do nb=0,8
             auxu   = auxu   + anb(nb,p)%intu
             auxv   = auxv   + anb(nb,p)%intv
             auxw   = auxw   + anb(nb,p)%intw
             auxomx = auxomx + anb(nb,p)%intomx
             auxomy = auxomy + anb(nb,p)%intomy
             auxomz = auxomz + anb(nb,p)%intomz
          enddo
          ap(p)%fxltot    = auxu
          ap(p)%fyltot    = auxv
          ap(p)%fzltot    = auxw
          ap(p)%torqxltot = auxomx
          ap(p)%torqyltot = auxomy
          ap(p)%torqzltot = auxomz
       endif
    enddo
    !$omp end parallel
    !
    return
  end subroutine updtlagrforces
  !
  subroutine updtintermediatevel(istep)
    implicit none
    integer,intent(in) :: istep
    integer i,j,k,p
    real sumfx,sumfy,sumfz
    real sumfx_all,sumfy_all,sumfz_all
    real forcextot_all,forceytot_all,forceztot_all
    real, parameter :: bulk_v_sup = 1.0  ! WARNING THIS NO LONGER WORKS
    !
    ! the spatial average of the particle-induced forces is subtracted, since
    ! we have periodic b.c.'s in all directions.
    ! (see Hoefler & Schwarzer, Phys. Rev. E, 61(6), 2000)
    !
    sumfx = 0.
    sumfy = 0.
    sumfz = 0.

    !$omp workshare
    sumfx = sum(forcex(1:imax,1:jmax,1:kmax))
    sumfy = sum(forcey(1:imax,1:jmax,1:kmax))
    sumfz = sum(forcez(1:imax,1:jmax,1:kmax))
    !$omp end workshare
    !!do k=1,kmax
    !!  do j=1,jmax
    !!    do i=1,imax
    !!      sumfx = sumfx + forcex(i,j,k)
    !!      sumfy = sumfy + forcey(i,j,k)
    !!      sumfz = sumfz + forcez(i,j,k)
    !!    enddo
    !!  enddo
    !!enddo

    call mpi_allreduce(sumfx,sumfx_all,1,mpi_real8,mpi_sum,comm_cart,error)
    call mpi_allreduce(sumfy,sumfy_all,1,mpi_real8,mpi_sum,comm_cart,error)
    call mpi_allreduce(sumfz,sumfz_all,1,mpi_real8,mpi_sum,comm_cart,error)
    sumfx = sumfx_all/(1.*itot*jtot*kmax)
    sumfy = sumfy_all/(1.*itot*jtot*kmax)
    sumfz = sumfz_all/(1.*itot*jtot*kmax)
    if (myid .eq. 0) then
       open(22,file='meanforces',position='append')
       write(22,'(I8,4E16.8)') istep,time/tref,sumfx,sumfy,sumfz
       close(22)
    endif
    !
    ! contribution from walls to forceytot
    !
    !
    ! Check: computation of averaged forces from lfps
    !
    forcextot = 0.
    forceytot = 0.
    forceztot = 0.
    !$omp parallel default(none) &
    !$omp& shared(ap,pmax) &
    !$omp& private(p)      &
    !$omp& reduction(+:forcextot,forceytot,forceztot)
    !$omp do
    do p=1,pmax
       if (ap(p)%mslv .gt. 0) then
          forcextot = forcextot + ap(p)%fxltot
          forceytot = forceytot + ap(p)%fyltot
          forceztot = forceztot + ap(p)%fzltot
       endif
    enddo
    !$omp end parallel
    call mpi_allreduce(forcextot,forcextot_all,1,mpi_real8,mpi_sum,comm_cart,error)
    call mpi_allreduce(forceytot,forceytot_all,1,mpi_real8,mpi_sum,comm_cart,error)
    call mpi_allreduce(forceztot,forceztot_all,1,mpi_real8,mpi_sum,comm_cart,error)
    forcextot = forcextot_all ! forcextot_all/(dVeul*itot*jtot*kmax) !averaged force per unit volume
    forceytot = forceytot_all ! forceytot_all/(dVeul*itot*jtot*kmax) !averaged force per unit volume
    forceztot = forceztot_all ! forceztot_all/(dVeul*itot*jtot*kmax) !averaged force per unit volume
    !
    if (myid .eq. 0) then
       open(23,file='IBM_forcing.out',position='append')
       write(23,'(I8,4E16.8)') istep,time,forcextot,forceytot,forceztot
       close(23)
      !  if ( (abs(sumfx-forcextot) .gt. 1.e-11) .or. (abs(sumfy-forceytot) .gt. 1.e-11) &
      !       .or. (abs(sumfz-forceztot) .gt. 1.e-11) ) then
      !     write(6,*) 'Error in computation of averaged particle forces!'
      !     write(6,*) 'Error in Fx = ',sumfx-forcextot
      !     write(6,*) 'Error in Fy = ',sumfy-forceytot
      !     write(6,*) 'Error in Fz = ',sumfz-forceztot
      !     !    call mpi_finalize(error)
      !     !    stop
      !  endif
    endif
   !  !
   !  ! flux imposed
   !  !
   !  !
   !  ! Forcing:
   !  !
   !  !forceytot = sumfy + wallshearnew
   !  !
   !  ! flux imposed
   !  !
   !  forcextot = 0.
   !  forceytot = 0.
   !  forceztot = 0.
   !  if(iniu.eq.'poi'.or.iniu.eq.'log') forceytot = -1.0*(bulk_v_sup-v_bulk)/dt
   !  !$omp parallel default(none) &
   !  !$omp&shared(dudt,dvdt,dwdt) &
   !  !$omp&shared(dudtold,dvdtold,dwdtold) &
   !  !$omp&shared(forcex,forcey,forcez) &
   !  !$omp&shared(forceytot,dt) &
   !  !$omp& private(i,j,k) 
   !  !$omp do
   !  do k=0,k1
   !    do j=0,j1
   !      do i=0,i1
   !        dudt(i,j,k) = dudtold(i,j,k) + forcex(i,j,k)*dt
   !        ! forceytot subtracted to get net zero forcing
   !        ! -dp/dy needed to balance total drag force = -forceytot
   !        dvdt(i,j,k) = dvdtold(i,j,k) + (forcey(i,j,k)-forceytot)*dt
   !        dwdt(i,j,k) = dwdtold(i,j,k) + forcez(i,j,k)*dt
   !      enddo
   !    enddo
   !  enddo
   !  !$omp end parallel
    !
  end subroutine updtintermediatevel
  !
end module mod_forcing
