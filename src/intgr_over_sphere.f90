module mod_intgr_over_sphere
  use mod_common
  use mod_common_mpi
  use mod_param
  use mod_tests
  implicit none
  private
  public intgr_over_sphere
contains
  !
!  subroutine intgr_over_sphere(intu,intv,intw,type) ! generalize for the three cases
  subroutine intgr_over_sphere(type)
    implicit none
    !real, intent(out), dimension(1:npmax) :: intu,intv,intw
    integer, intent(in) :: type
    type pneighbor
       real :: x,y,z, &
            intu,intv,intw
    end type pneighbor
    type(pneighbor), dimension(0:8,1:npmax) :: anb ! can be pmax because it is a subroutine-specific array!
    
    integer :: i,j,k,p,cp
    integer :: ilow,ihigh,jlow,jhigh,klow,khigh
    real :: boundleftnb,boundrightnb,boundfrontnb,boundbacknb
    real :: coorxc,cooryc,coorzc
    real :: coorxmin,coorxplus,coorymin,cooryplus,coorzmin,coorzplus
    real :: radin,radin2,dist2
    real :: dx,dy,dz
    real :: sum1,sum2
    integer :: nb,nbsend,nbrecv
    integer :: nrrequests ! MPI request counter (SS)
    integer :: arrayrequests(1:3) ! 3=3*1 (master might have 3 slaves)
    integer :: arraystatuses(MPI_STATUS_SIZE,1:3)
    real  :: sgndist(1:8)
    real :: rx,ry,rz
    integer :: idp,tag
    character(len=5) rankpr
    real :: dxm2,dxp2,dym2,dyp2,dzm2,dzp2
    real :: auxu,auxv,auxw
    !
    dx       = 1./dxi
    dy       = 1./dyi
    dz       = 1./dzi
    radin    = radius-2.*dx
    radin2   = radin**2

    if (runmode == RUNMODE_TEST_PAR_FREE_FALL_NO_SCALAR) then
       call test_sphere_int_param(radin, radin2, dx, dy, dz,type)
    end if
    
    !
    ! for r < radius - sqrt( (dx)**2 + (dy)**2 + (dz)**2 ) = 
    !        radius - 1.732050808*dx : all cell corner points within sphere
    !
    ! first step: slaves need from master x,..,xfp,.. etc
    ! anb(1:8,p) will contain x,y,z of particle p as seen from neighbor 1..8 (SS)
    !$omp workshare
    anb(0,1:npmax)%x = ap(1:npmax)%x ! ap: 'a particle' array, 
    anb(0,1:npmax)%y = ap(1:npmax)%y
    anb(0,1:npmax)%z = ap(1:npmax)%z


    forall(p=1:npmax)
       anb(1:8,p)%x = 0.
       anb(1:8,p)%y = 0.
       anb(1:8,p)%z = 0.
    end forall
    !$omp end workshare
    ! 
    do p=1,pmax
       nrrequests = 0  ! Initialize MPI request counter
       do nb=1,8
          nbsend = nb    ! neighbor(nbsend) = rank of process to which data is send
          idp = abs(ap(p)%mslv)
          tag = idp*10+nbsend-idp*10
          nbrecv  = nb+4  ! neighbor(nbrecv) = rank of process from which data is received
          if (nbrecv .gt. 8) nbrecv = nbrecv - 8
          if (ap(p)%mslv .gt. 0) then 
             ! myid is master of particle ap(p)%mslv
             if (ap(p)%nb(nbsend) .eq. 1) then
                if (neighbor(nbsend).eq.myid) then
                   ! process might be both master and slave of same particle due to periodic b.c.'s
                   anb(nbrecv,p)%x = anb(0,p)%x
                   anb(nbrecv,p)%y = anb(0,p)%y
                   anb(nbrecv,p)%z = anb(0,p)%z
                   ! recompute particle positions due to periodic b.c.'s
                   if (nbrecv .eq. 1) then
                      anb(nbrecv,p)%x = anb(nbrecv,p)%x + lx
                   endif
                   if (nbrecv .eq. 2) then
                      anb(nbrecv,p)%x = anb(nbrecv,p)%x + lx
                      anb(nbrecv,p)%y = anb(nbrecv,p)%y - ly
                   endif
                   if (nbrecv .eq. 3) then
                      anb(nbrecv,p)%y = anb(nbrecv,p)%y - ly
                   endif
                   if (nbrecv .eq. 4) then
                      anb(nbrecv,p)%x = anb(nbrecv,p)%x - lx
                      anb(nbrecv,p)%y = anb(nbrecv,p)%y - ly
                   endif
                   if (nbrecv .eq. 5) then
                      anb(nbrecv,p)%x = anb(nbrecv,p)%x - lx
                   endif
                   if (nbrecv .eq. 6) then
                      anb(nbrecv,p)%x = anb(nbrecv,p)%x - lx
                      anb(nbrecv,p)%y = anb(nbrecv,p)%y + ly
                   endif
                   if (nbrecv .eq. 7) then
                      anb(nbrecv,p)%y = anb(nbrecv,p)%y + ly
                   endif
                   if (nbrecv .eq. 8) then
                      anb(nbrecv,p)%x = anb(nbrecv,p)%x + lx
                      anb(nbrecv,p)%y = anb(nbrecv,p)%y + ly
                   endif
                else
                   ! neighbor(nbsend) is rank of slave for particle ap(p)%mslv
                   nrrequests = nrrequests + 1
                   call MPI_ISEND(anb(0,p)%x,3,MPI_REAL8,neighbor(nbsend), &
                        tag,comm_cart,arrayrequests((nrrequests-1) + 1),error)
                   ! send anb(0,p) x,y,z data to neighbor(nbsend), use tag to distinguish (SS)
                   ! send x,y,z -> 3 contiguous info
                   ! (see definition of type pneighbor in the begining of the subroutine)
                endif
             endif
          endif
          if (ap(p)%mslv .lt. 0) then
             ! myid is slave of particle -ap(p)%mslv
             if (ap(p)%nb(nbrecv) .eq. 1) then
                ! neighbor(nbrecv) is rank of master of particle -ap(p)%mslv
                nrrequests = nrrequests + 1
                call MPI_IRECV(anb(nbrecv,p)%x,3,MPI_REAL8,neighbor(nbrecv), &
                     tag,comm_cart,arrayrequests((nrrequests-1) + 1),error)
                ! recv anb(nbrecv,p) x,y,z data from neighbor(nbrecv), use tag to distinguish (SS)
                ! recv x,y,z -> 3 contiguous info
                ! (see definition of type pneighbor in the begining of the subroutine)
             endif
          endif
       enddo ! do nb=
       call MPI_WAITALL(nrrequests,arrayrequests,arraystatuses,error)
    enddo
    
   !  ! Debug sphere integration--check first step results
   !  ! Have to put here because the anb have three different definitions
   !  ! in forcing.f90, phase_indicator.f90 and intgr_over_sphere.f90.    
   !  ! for debugging sphere integration issues (SS)
   ! if (runmode == RUNMODE_TEST_PAR_FREE_FALL_NO_SCALAR) then     
   !  if (myid == 6 .and. type == 1) then
   !     write(6,*) ''
   !     write(6,*) '=== INT SPHERE FIRST STEP TEST ==='
   !     do p = 1, min(pmax, 3)  ! Check first few particles
   !        if (ap(p)%mslv .ne. 0) then
   !           write(6,*) '--- Particle', p, '---'
   !           write(6,'(A,I0)') 'myid = ', myid
   !           write(6,'(A,I0,A,3F12.6)') 'ap(',p,')%x,y,z = ', ap(p)%x, ap(p)%y, ap(p)%z
   !           write(6,'(A,I0,A,3F12.6)') 'anb(0,',p,')%x,y,z = ', anb(0,p)%x, anb(0,p)%y, anb(0,p)%z
   !           write(6,'(A,I0,A,I0)') 'ap(',p,')%mslv = ', ap(p)%mslv
   !           do nb = 1, 8
   !                 write(6,'(A,I0,A,I0,A,3F12.6)') 'anb(',nb,',',p,')%x,y,z = ', anb(nb,p)%x, anb(nb,p)%y, anb(nb,p)%z
   !           enddo
   !        endif
   !     enddo
   !     write(6,*) '=== END FIRST STEP TEST ==='
   !     write(6,*) ''
   !  end if    
   ! end if
    !
    ! second step: recompute particle positions for slaves due to periodic b.c.'s.
    ! Required: (part of) particle within domain bounds of slave process.
    !
    !$omp parallel default(none) &
    !$omp&shared(ap,anb,pmax,coords)  &
    !$omp&private(p,nbrecv,boundleftnb,boundbacknb,boundrightnb,boundfrontnb)  
    !$omp do
    do p=1,pmax
       if (ap(p)%mslv .lt. 0) then
          ! myid is slave of particle -ap(p)%mslv
          nbrecv=1  ! Taking the slave process as reference (SS)
          if (ap(p)%nb(nbrecv) .eq. 1) then
             ! neighbor(nbrecv) is rank of master of particle -ap(p)%mslv
             boundleftnb  = (coords(1)+1)*lx/(1.*dims(1)) ! left boundary of neighbor nb
             if (anb(nbrecv,p)%x .lt. boundleftnb) then
                anb(nbrecv,p)%x = anb(nbrecv,p)%x + lx
             endif
          endif
          nbrecv=2
          if (ap(p)%nb(nbrecv) .gt. 0) then
             boundleftnb  = (coords(1)+1)*lx/(1.*dims(1)) ! left boundary of neighbor nb
             boundbacknb  = (coords(2))*ly/(1.*dims(2)) ! back boundary of neighbor nb
             if (anb(nbrecv,p)%x .lt. boundleftnb) then
                anb(nbrecv,p)%x = anb(nbrecv,p)%x + lx
             endif
             if (anb(nbrecv,p)%y .gt. boundbacknb) then
                anb(nbrecv,p)%y = anb(nbrecv,p)%y - ly
             endif
          endif
          nbrecv=3
          if (ap(p)%nb(nbrecv) .gt. 0) then
             boundbacknb  = (coords(2))*ly/(1.*dims(2)) ! back boundary of neighbor nb
             if (anb(nbrecv,p)%y .gt. boundbacknb) then
                anb(nbrecv,p)%y = anb(nbrecv,p)%y - ly
             endif
          endif
          nbrecv=4
          if (ap(p)%nb(nbrecv) .gt. 0) then
             boundrightnb = (coords(1))*lx/(1.*dims(1)) ! right boundary of neighbor nb
             boundbacknb  = (coords(2))*ly/(1.*dims(2)) ! back  boundary of neighbor nb
             if (anb(nbrecv,p)%x .gt. boundrightnb) then
                anb(nbrecv,p)%x = anb(nbrecv,p)%x - lx
             endif
             if (anb(nbrecv,p)%y .gt. boundbacknb) then
                anb(nbrecv,p)%y = anb(nbrecv,p)%y - ly
             endif
          endif
          nbrecv=5
          if (ap(p)%nb(nbrecv) .gt. 0) then
             boundrightnb = (coords(1))*lx/(1.*dims(1)) ! right boundary of neighbor nb
             if (anb(nbrecv,p)%x .gt. boundrightnb) then
                anb(nbrecv,p)%x = anb(nbrecv,p)%x - lx
             endif
          endif
          nbrecv=6
          if (ap(p)%nb(nbrecv) .gt. 0) then
             boundrightnb = (coords(1))*lx/(1.*dims(1)) ! right boundary of neighbor nb
             boundfrontnb = (coords(2)+1)*ly/(1.*dims(2)) ! front boundary of neighbor nb
             if (anb(nbrecv,p)%x .gt. boundrightnb) then
                anb(nbrecv,p)%x = anb(nbrecv,p)%x - lx
             endif
             if (anb(nbrecv,p)%y .lt. boundfrontnb) then
                anb(nbrecv,p)%y = anb(nbrecv,p)%y + ly
             endif
          endif
          nbrecv=7
          if (ap(p)%nb(nbrecv) .gt. 0) then
             boundfrontnb = (coords(2)+1)*ly/(1.*dims(2)) ! front boundary of neighbor nb
             if (anb(nbrecv,p)%y .lt. boundfrontnb) then
                anb(nbrecv,p)%y = anb(nbrecv,p)%y + ly
             endif
          endif
          nbrecv=8
          if (ap(p)%nb(nbrecv) .gt. 0) then
             boundleftnb  = (coords(1)+1)*lx/(1.*dims(1)) ! left  boundary of neighbor nb
             boundfrontnb = (coords(2)+1)*ly/(1.*dims(2)) ! front boundary of neighbor nb
             if (anb(nbrecv,p)%x .lt. boundleftnb) then
                anb(nbrecv,p)%x = anb(nbrecv,p)%x + lx
             endif
             if (anb(nbrecv,p)%y .lt. boundfrontnb) then
                anb(nbrecv,p)%y = anb(nbrecv,p)%y + ly
             endif
          endif
       endif
    enddo
    !$omp end parallel
    !
   !  ! Debug sphere integration--check first step results
   !  ! Have to put here because the anb have three different definitions
   !  ! in forcing.f90, phase_indicator.f90 and intgr_over_sphere.f90.    
   !  ! for debugging sphere integration issues (SS)
   ! if (runmode == RUNMODE_TEST_PAR_FREE_FALL_NO_SCALAR) then     
   !  if (myid == 9 .and. type == 1) then
   !     write(6,*) ''
   !     write(6,*) '=== INT SPHERE SECOND STEP PART 1 TEST ==='
   !     do p = 1, min(pmax, 3)  ! Check first few particles
   !        if (ap(p)%mslv .ne. 0) then
   !           write(6,*) '--- Particle', p, '---'
   !           write(6,'(A,I0)') 'myid = ', myid
   !           write(6,'(A,I0,A,3F12.6)') 'ap(',p,')%x,y,z = ', ap(p)%x, ap(p)%y, ap(p)%z
   !           write(6,'(A,I0,A,3F12.6)') 'anb(0,',p,')%x,y,z = ', anb(0,p)%x, anb(0,p)%y, anb(0,p)%z
   !           write(6,'(A,I0,A,I0)') 'ap(',p,')%mslv = ', ap(p)%mslv
   !           do nb = 1, 8
   !                 write(6,'(A,I0,A,I0,A,3F12.6)') 'anb(',nb,',',p,')%x,y,z = ', anb(nb,p)%x, anb(nb,p)%y, anb(nb,p)%z
   !           enddo
   !        endif
   !     enddo
   !     write(6,*) '=== END SECOND STEP PART 1 TEST ==='
   !     write(6,*) ''
   !  end if   
   ! end if
    ! second step: perform integration.
    !
    !$omp workshare
    forall(p=1:npmax,i=0:8) 
       anb(i,p)%intu = 0.
       anb(i,p)%intv = 0.
       anb(i,p)%intw = 0.
    end forall
    !$omp end workshare
    !
    !$omp parallel default(none) &
    !$omp&shared(ap,anb,pmax) & 
    !$omp&shared(unew,vnew,wnew) &
    !$omp&shared(type,dx,dy,dz,radin2,dveul) &
    !$omp&shared(myid,neighbor,boundleftmyid,boundfrontmyid) &
    !$omp&private(p,nb,nbrecv,coorxc,cooryc,coorzc,ilow,ihigh,jlow,jhigh,klow,khigh) &
    !$omp&private(coorxmin,coorxplus,coorymin,cooryplus,coorzmin,coorzplus) &
    !$omp&private(dxm2,dxp2,dym2,dyp2,dzm2,dzp2,dist2,rx,ry,rz,sgndist,sum1,sum2) &
    !$omp&private(i,j,k) 
    !$omp do
    

    do p=1,pmax
       if (ap(p)%mslv .ne. 0) then
          ! myid is master or slave of particle abs(ap(p)%mslv)
          do nb=0,8
             if(( (nb.gt.0) .and. (ap(p)%mslv.lt.0) .and. (ap(p)%nb(nb).eq.1) ) .or. & !slave
                ( (nb.eq.0) .and. (ap(p)%mslv.gt.0) ) .or. & !pure master
                ( (nb.gt.0) .and. (ap(p)%mslv.gt.0) .and. (ap(p)%nb(nb).eq.1) .and. (neighbor(nb).eq.myid) ) ) then
                ! master that looks like a slave due to periodic bcs
                nbrecv = nb
                ! neighbor(nbrecv) is rank of master of particle -ap(p)%mslv
                if ( neighbor(nb) .eq. myid .and. nb.gt.0.and.ap(p)%mslv.gt.0) then
                   nbrecv = nb + 4
                   if (nbrecv .gt. 8) nbrecv = nbrecv-8
                endif
                coorxc = anb(nbrecv,p)%x-boundleftmyid
                cooryc = anb(nbrecv,p)%y-boundfrontmyid
                coorzc = anb(nbrecv,p)%z
                ilow  = nint( (coorxc-radius)*dxi - 2. )
                ihigh = nint( (coorxc+radius)*dxi + 2. )
                jlow  = nint( (cooryc-radius)*dyi - 2. )
                jhigh = nint( (cooryc+radius)*dyi + 2. )
                klow  = nint( (coorzc-radius)*dzi - 2. )
                khigh = nint( (coorzc+radius)*dzi + 2. )
                if (runmode == RUNMODE_TEST_PAR_FREE_FALL_NO_SCALAR) then
                   call test_sphere_int_second_part2(nb, p, coorxc,cooryc,coorzc,ilow,ihigh,jlow,jhigh,klow,khigh,type)
                end if  
                if (ilow .lt. 1) ilow = 1
                if (jlow .lt. 1) jlow = 1
                if (klow .lt. 1) klow = 1
                if (ihigh .gt. imax) ihigh = imax
                if (jhigh .gt. jmax) jhigh = jmax
                if (khigh .gt. kmax) khigh = kmax
                !
                ! u-velocity
                !
                do k=klow,khigh
                   coorzmin  = (k-1)*dz
                   coorzplus = k*dz
                   dzp2 = (coorzplus - anb(nbrecv,p)%z)**2.
                   do j=jlow,jhigh
                      coorymin  = boundfrontmyid + (j-1)*dy
                      cooryplus = boundfrontmyid + j*dy
                      dyp2 = (cooryplus - anb(nbrecv,p)%y)**2.
                      do i=ilow,ihigh
                         coorxmin  = boundleftmyid + (i-0.5)*dx
                         coorxplus = boundleftmyid + (i+0.5)*dx
                         dxp2 = (coorxplus - anb(nbrecv,p)%x)**2.
                         dist2 = dxp2+dyp2+dzp2 
                         if (dist2 .lt. radin2) then
                            if (type.eq.1) then
                               anb(nbrecv,p)%intu = anb(nbrecv,p)%intu + dVeul
                            elseif(type.eq.2) then
                               anb(nbrecv,p)%intu = anb(nbrecv,p)%intu + unew(i,j,k)*dVeul
                            elseif(type.eq.3) then
                               ry = boundfrontmyid + (j-0.5)*dy - anb(nbrecv,p)%y
                               rz = (k-0.5)*dz                  - anb(nbrecv,p)%z
                               anb(nbrecv,p)%intv = anb(nbrecv,p)%intv + rz*unew(i,j,k)*dVeul
                               anb(nbrecv,p)%intw = anb(nbrecv,p)%intw - ry*unew(i,j,k)*dVeul
                            endif
                         else
                            dxm2 = (coorxmin - anb(nbrecv,p)%x)**2.
                            dym2 = (coorymin - anb(nbrecv,p)%y)**2.
                            dzm2 = (coorzmin - anb(nbrecv,p)%z)**2.
                            sgndist(1) = radius - sqrt(dxm2+dym2+dzm2) ! left-front-bottom corner
                            sgndist(2) = radius - sqrt(dxm2+dyp2+dzm2) ! left-back-bottom corner
                            sgndist(3) = radius - sqrt(dxp2+dyp2+dzm2) ! right-back-bottom corner
                            sgndist(4) = radius - sqrt(dxp2+dym2+dzm2) ! right-front-bottom corner
                            sgndist(5) = radius - sqrt(dxm2+dym2+dzp2) ! left-front-top corner
                            sgndist(6) = radius - sqrt(dxm2+dyp2+dzp2) ! left-back-top corner
                            sgndist(7) = radius - sqrt(dxp2+dyp2+dzp2) ! right-back-top corner
                            sgndist(8) = radius - sqrt(dxp2+dym2+dzp2) ! right-front-top corner
                            sum1 = 0.
                            sum2 = 0.
                            do cp=1,8 ! cp = corner point
                               if (sgndist(cp) .gt. 0.) sum1 = sum1 + sgndist(cp)
                               sum2 = sum2 + abs( sgndist(cp) )
                            enddo
                            if (type.eq.1) then
                               anb(nbrecv,p)%intu = anb(nbrecv,p)%intu + (sum1/sum2)*dVeul
                            elseif(type.eq.2) then
                               anb(nbrecv,p)%intu = anb(nbrecv,p)%intu + unew(i,j,k)*(sum1/sum2)*dVeul
                            elseif(type.eq.3) then
                               ry = boundfrontmyid + (j-0.5)*dy - anb(nbrecv,p)%y
                               rz = (k-0.5)*dz                  - anb(nbrecv,p)%z
                               anb(nbrecv,p)%intv = anb(nbrecv,p)%intv + rz*unew(i,j,k)*(sum1/sum2)*dVeul
                               anb(nbrecv,p)%intw = anb(nbrecv,p)%intw - ry*unew(i,j,k)*(sum1/sum2)*dVeul
                            endif
                         endif
                      enddo ! do i=
                   enddo ! do j=
                enddo ! do k=
                !
                ! v-velocity
                !
                do k=klow,khigh
                   coorzmin  = (k-1)*dz
                   coorzplus = k*dz
                   dzp2 = (coorzplus - anb(nbrecv,p)%z)**2.
                   do j=jlow,jhigh
                      coorymin  = boundfrontmyid + (j-0.5)*dy
                      cooryplus = boundfrontmyid + (j+0.5)*dy
                      dyp2 = (cooryplus - anb(nbrecv,p)%y)**2.
                      do i=ilow,ihigh
                         coorxmin  = boundleftmyid + (i-1)*dx
                         coorxplus = boundleftmyid + i*dx
                         dxp2 = (coorxplus - anb(nbrecv,p)%x)**2.
                         dist2 = dxp2+dyp2+dzp2 
                         if (dist2 .lt. radin2) then
                            if (type.eq.1) then
                               anb(nbrecv,p)%intv = anb(nbrecv,p)%intv + dVeul
                            elseif(type.eq.2) then
                               anb(nbrecv,p)%intv = anb(nbrecv,p)%intv + vnew(i,j,k)*dVeul
                            elseif(type.eq.3) then
                               rx = boundleftmyid + (i-0.5)*dx - anb(nbrecv,p)%x
                               rz = (k-0.5)*dz                 - anb(nbrecv,p)%z
                               anb(nbrecv,p)%intu = anb(nbrecv,p)%intu - rz*vnew(i,j,k)*dVeul
                               anb(nbrecv,p)%intw = anb(nbrecv,p)%intw + rx*vnew(i,j,k)*dVeul
                            endif
                         else
                            dxm2 = (coorxmin - anb(nbrecv,p)%x)**2.
                            dym2 = (coorymin - anb(nbrecv,p)%y)**2.
                            dzm2 = (coorzmin - anb(nbrecv,p)%z)**2.
                            sgndist(1) = radius - sqrt(dxm2+dym2+dzm2) ! left-front-bottom corner
                            sgndist(2) = radius - sqrt(dxm2+dyp2+dzm2) ! left-back-bottom corner
                            sgndist(3) = radius - sqrt(dxp2+dyp2+dzm2) ! right-back-bottom corner
                            sgndist(4) = radius - sqrt(dxp2+dym2+dzm2) ! right-front-bottom corner
                            sgndist(5) = radius - sqrt(dxm2+dym2+dzp2) ! left-front-top corner
                            sgndist(6) = radius - sqrt(dxm2+dyp2+dzp2) ! left-back-top corner
                            sgndist(7) = radius - sqrt(dxp2+dyp2+dzp2) ! right-back-top corner
                            sgndist(8) = radius - sqrt(dxp2+dym2+dzp2) ! right-front-top corner
                            sum1 = 0.
                            sum2 = 0.
                            do cp=1,8 ! cp = corner point
                               if (sgndist(cp) .gt. 0.) sum1 = sum1 + sgndist(cp)
                               sum2 = sum2 + abs( sgndist(cp) )
                            enddo
                            if (type.eq.1) then
                               anb(nbrecv,p)%intv = anb(nbrecv,p)%intv + (sum1/sum2)*dVeul
                            elseif(type.eq.2) then
                               anb(nbrecv,p)%intv = anb(nbrecv,p)%intv + vnew(i,j,k)*(sum1/sum2)*dVeul
                            elseif(type.eq.3) then
                               rx = boundleftmyid + (i-0.5)*dx - anb(nbrecv,p)%x
                               rz = (k-0.5)*dz                 - anb(nbrecv,p)%z 
                               anb(nbrecv,p)%intu = anb(nbrecv,p)%intu - rz*vnew(i,j,k)*(sum1/sum2)*dVeul
                               anb(nbrecv,p)%intw = anb(nbrecv,p)%intw + rx*vnew(i,j,k)*(sum1/sum2)*dVeul
                            endif
                         endif
                      enddo ! do i=
                   enddo ! do j=
                enddo ! do k=
                !
                ! w-velocity
                !
                do k=klow,khigh
                   coorzmin  = (k-0.5)*dz
                   coorzplus = (k+0.5)*dz
                   dzp2 = (coorzplus - anb(nbrecv,p)%z)**2.
                   do j=jlow,jhigh
                      coorymin  = boundfrontmyid + (j-1)*dy
                      cooryplus = boundfrontmyid + j*dy
                      dyp2 = (cooryplus - anb(nbrecv,p)%y)**2.
                      do i=ilow,ihigh
                         coorxmin  = boundleftmyid + (i-1)*dx
                         coorxplus = boundleftmyid + i*dx
                         dxp2 = (coorxplus - anb(nbrecv,p)%x)**2.
                         dist2 = dxp2+dyp2+dzp2 
                         if (dist2 .lt. radin2) then
                            if (type.eq.1) then
                               anb(nbrecv,p)%intw = anb(nbrecv,p)%intw + dVeul
                            elseif(type.eq.2) then
                               anb(nbrecv,p)%intw = anb(nbrecv,p)%intw + wnew(i,j,k)*dVeul
                            elseif(type.eq.3) then
                               rx = boundleftmyid  + (i-0.5)*dx - anb(nbrecv,p)%x
                               ry = boundfrontmyid + (j-0.5)*dy - anb(nbrecv,p)%y
                               anb(nbrecv,p)%intu = anb(nbrecv,p)%intu + ry*wnew(i,j,k)*dVeul
                               anb(nbrecv,p)%intv = anb(nbrecv,p)%intv - rx*wnew(i,j,k)*dVeul
                            endif
                         else
                            dxm2 = (coorxmin - anb(nbrecv,p)%x)**2.
                            dym2 = (coorymin - anb(nbrecv,p)%y)**2.
                            dzm2 = (coorzmin - anb(nbrecv,p)%z)**2.
                            sgndist(1) = radius - sqrt(dxm2+dym2+dzm2) ! left-front-bottom corner
                            sgndist(2) = radius - sqrt(dxm2+dyp2+dzm2) ! left-back-bottom corner
                            sgndist(3) = radius - sqrt(dxp2+dyp2+dzm2) ! right-back-bottom corner
                            sgndist(4) = radius - sqrt(dxp2+dym2+dzm2) ! right-front-bottom corner
                            sgndist(5) = radius - sqrt(dxm2+dym2+dzp2) ! left-front-top corner
                            sgndist(6) = radius - sqrt(dxm2+dyp2+dzp2) ! left-back-top corner
                            sgndist(7) = radius - sqrt(dxp2+dyp2+dzp2) ! right-back-top corner
                            sgndist(8) = radius - sqrt(dxp2+dym2+dzp2) ! right-front-top corner
                            sum1 = 0.
                            sum2 = 0.
                            do cp=1,8 ! cp = corner point
                               if (sgndist(cp) .gt. 0.) sum1 = sum1 + sgndist(cp)
                               sum2 = sum2 + abs( sgndist(cp) )
                            enddo
                            if (type.eq.1) then
                               anb(nbrecv,p)%intw = anb(nbrecv,p)%intw + (sum1/sum2)*dVeul
                            elseif(type.eq.2) then
                               anb(nbrecv,p)%intw = anb(nbrecv,p)%intw + wnew(i,j,k)*(sum1/sum2)*dVeul
                            elseif(type.eq.3) then
                               rx = boundleftmyid  + (i-0.5)*dx - anb(nbrecv,p)%x
                               ry = boundfrontmyid + (j-0.5)*dy - anb(nbrecv,p)%y 
                               anb(nbrecv,p)%intu = anb(nbrecv,p)%intu + ry*wnew(i,j,k)*(sum1/sum2)*dVeul
                               anb(nbrecv,p)%intv = anb(nbrecv,p)%intv - rx*wnew(i,j,k)*(sum1/sum2)*dVeul
                            endif
                         endif
                      enddo ! do i=
                   enddo ! do j=
                enddo ! do k=
             endif
          enddo ! do nbrecv=
       endif
    enddo
    !$omp end parallel
    !

    ! third step: communicate data of slaves to their masters
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
                   call MPI_IRECV(anb(nbsend,p)%intu,3,MPI_REAL8,neighbor(nbsend), &
                        tag,comm_cart,arrayrequests((nrrequests-1) + 1),error)
                   ! recv intu,intv,intwx,y,z -> 3 contiguous info
                   ! (see definition of type pneighbor in the begining of the subroutine)
                endif
             endif
          endif
          if (ap(p)%mslv .lt. 0) then
             ! myid is slave of particle -ap(p)%mslv
             if (ap(p)%nb(nbrecv) .eq. 1) then
                ! neighbor(nbrecv) is rank of master of particle -ap(p)%mslv
                nrrequests = nrrequests + 1
                call MPI_ISEND(anb(nbrecv,p)%intu,3,MPI_REAL8,neighbor(nbrecv), &
                     tag,comm_cart,arrayrequests((nrrequests-1) + 1),error)
                ! send intu,intv,intwx,y,z -> 3 contiguous info
                ! (see definition of type pneighbor in the begining of the subroutine)
             endif
          endif
       enddo ! do nb=
       call MPI_WAITALL(nrrequests,arrayrequests,arraystatuses,error)
    enddo
    !
    ! Sum all contributions together.
    !
    !$omp parallel default(none) &
    !$omp&shared(ap,anb,pmax,type) &
    !$omp&private(p,nb) reduction(+:auxu,auxv,auxw)
    SELECT CASE (type)
    CASE (1,2)
       !$omp do 
       do p=1,pmax
          if (ap(p)%mslv .gt. 0) then
             auxu = 0.
             auxv = 0.
             auxw = 0.
             do nb=0,8
                auxu = auxu + anb(nb,p)%intu
                auxv = auxv + anb(nb,p)%intv
                auxw = auxw + anb(nb,p)%intw
             enddo
             ap(p)%intu = auxu
             ap(p)%intv = auxv
             ap(p)%intw = auxw
          else
             ap(p)%intu = 0.
             ap(p)%intv = 0.
             ap(p)%intw = 0.
          end if
       enddo
    CASE (3)
       !$omp do 
       do p=1,pmax
          if (ap(p)%mslv .gt. 0) then
             auxu = 0.
             auxv = 0.
             auxw = 0.
             do nb=0,8
                auxu = auxu + anb(nb,p)%intu
                auxv = auxv + anb(nb,p)%intv
                auxw = auxw + anb(nb,p)%intw
             enddo
             ap(p)%intomx = auxu
             ap(p)%intomy = auxv
             ap(p)%intomz = auxw
          else
             ap(p)%intomx = 0.
             ap(p)%intomy = 0.
             ap(p)%intomz = 0.
          end if
       enddo
    END SELECT
    !$omp end parallel
    !

    if (type.eq.1) then
       write(rankpr,'(i5.5)') myid
       do p=1,pmax
          if (ap(p)%mslv .gt. 0) then
             idp = ap(p)%mslv
             !    open(22,file=datadir//'volsphr'//rankpr//'.txt',position='append')
             !    write(22,'(2I8,4E16.8)') myid,idp,intu(p),intv(p),intw(p),Volp
             !    close(22)
             if (idp .eq. 1) then
                write(6,'(A31,E16.8)') 'Calc. value of volume sphere = ',ap(p)%intu
                write(6,'(A31,E16.8)') 'Exact value of volume sphere = ',Volp
                write(6,'(A41,E16.8)') 'Error in calc. of volume sphere (in %) = ',100.*(ap(p)%intu-Volp)/Volp
             endif
          endif
       enddo
    endif
    !
    return
  end subroutine intgr_over_sphere
  !
end module mod_intgr_over_sphere
