module mod_phase_indicator
  use mod_common
  use mod_common_mpi
  use mod_param
  use mod_bound
  implicit none
  private
  public phase_indicator
contains
  !
  subroutine phase_indicator(gamu,gamv,gamw,gamp, &
       upart,vpart,wpart,itype)
    implicit none
    real, intent(out), dimension(0:i1,0:j1,0:k1) :: gamu,gamv,gamw,gamp,upart,vpart,wpart 
    integer, intent(in) :: itype ! 1 -> impose rigid body motion -> translation + rotation
    ! 0 -> impose rigid body motion -> rotation
    type pneighbor
       real :: x,y,z,u,v,w,omx,omy,omz
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
    integer :: nrrequests
    integer :: arrayrequests(1:3) ! 3=3*1 (master might have 3 slaves)
    integer :: arraystatuses(MPI_STATUS_SIZE,1:3)
    real  :: sgndist(1:8)
    real :: rx,ry,rz
    integer :: idp,tag
    character(len=3) rankpr
    real :: dxm2,dxp2,dym2,dyp2,dzm2,dzp2
    real :: aux,aux_all 
    !
    dx       = 1./dxi
    dy       = 1./dyi
    dz       = 1./dzi
    radin    = radius-2.*dx
    radin2   = radin**2
    !
    ! for r < radius - sqrt( (dx)**2 + (dy)**2 + (dz)**2 ) = 
    !        radius - 1.732050808*dx : all cell corner points within sphere
    !
    ! first step: slaves need from master x,..,xfp,.. etc
    !
    !$omp workshare
    anb(0,1:pmax)%x = ap(1:pmax)%x
    anb(0,1:pmax)%y = ap(1:pmax)%y
    anb(0,1:pmax)%z = ap(1:pmax)%z
    anb(0,1:pmax)%u = ap(1:pmax)%u
    anb(0,1:pmax)%v = ap(1:pmax)%v
    anb(0,1:pmax)%w = ap(1:pmax)%w
    anb(0,1:pmax)%omx = ap(1:pmax)%omx
    anb(0,1:pmax)%omy = ap(1:pmax)%omy
    anb(0,1:pmax)%omz = ap(1:pmax)%omz
    forall(p=1:pmax)
       anb(1:8,p)%x = 0.
       anb(1:8,p)%y = 0.
       anb(1:8,p)%z = 0.
       anb(1:8,p)%u = 0.
       anb(1:8,p)%v = 0.
       anb(1:8,p)%w = 0.
       anb(1:8,p)%omx = 0.
       anb(1:8,p)%omy = 0.
       anb(1:8,p)%omz = 0.
    end forall
    !$omp end workshare
    !
    do p=1,pmax
       nrrequests = 0
       do nb=1,8
          nbsend = nb    ! neighbor(nbsend) = rank of process to which data is send
          idp = abs(ap(p)%mslv)
          tag = idp
          nbrecv  = nb+4  ! neighbor(nbrecv)  = rank of process from which data is received
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
                   call MPI_ISEND(anb(0,p)%x,9,MPI_REAL8,neighbor(nbsend), &
                        tag,comm_cart,arrayrequests((nrrequests-1) + 1),error)
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
                call MPI_IRECV(anb(nbrecv,p)%x,9,MPI_REAL8,neighbor(nbrecv), &
                     tag,comm_cart,arrayrequests((nrrequests-1) + 1),error)
                ! recv x,y,z -> 3 contiguous info
                ! (see definition of type pneighbor in the begining of the subroutine)
             endif
          endif
       enddo ! do nb=
       call MPI_WAITALL(nrrequests,arrayrequests,arraystatuses,error)
    enddo
    !
    ! second step: recompute particle positions for slaves due to periodic b.c.'s.
    ! Required: (part of) particle within domain bounds of slave process.
    !
    !$omp parallel default(shared) &
    !$omp&private(p,nbrecv,boundleftnb,boundbacknb,boundrightnb,boundfrontnb)  
    !$omp do
    do p=1,pmax
       if (ap(p)%mslv .lt. 0) then
          ! myid is slave of particle -ap(p)%mslv
          nbrecv=1
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
    ! third step: perform integration.
    !
    gamp(:,:,:) = 0.
    gamu(:,:,:) = 0.
    gamv(:,:,:) = 0.
    gamw(:,:,:) = 0.
    upart(:,:,:) = 0.
    vpart(:,:,:) = 0.
    wpart(:,:,:) = 0.
    !$omp parallel default(shared) &
    !$omp&private(p,nb,nbrecv,coorxc,cooryc,coorzc,ilow,ihigh,jlow,jhigh,klow,khigh) &
    !$omp&private(coorxmin,coorxplus,coorymin,cooryplus,coorzmin,coorzplus) &
    !$omp&private(dxm2,dxp2,dym2,dyp2,dzm2,dzp2,dist2,rx,ry,rz,sgndist,sum1,sum2) &
    !$omp&private(i,j,k) 
    !$omp do 
    do p=1,pmax
       if (ap(p)%mslv .ne. 0) then
          ! myid is master or slave of particle abs(ap(p)%mslv)
          do nb=0,8
             if((nb.gt.0.and.ap(p)%mslv.lt.0.and.ap(p)%nb(nb).eq.1) .or. & !slave
                  (nb.eq.0.and.ap(p)%mslv.gt.0) .or. & !pure master
                  (nb.gt.0.and.ap(p)%mslv.gt.0.and.ap(p)%nb(nb).eq.1.and.neighbor(nb) .eq. myid)) then
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
                         if (dist2 .lt. radin2.and.itype.eq.1) then
                            gamu(i,j,k) = 1. 
                            upart(i,j,k) = anb(nbrecv,p)%u + &
                                 anb(nbrecv,p)%omy*((k-0.5)*dz-anb(nbrecv,p)%z) - &
                                 anb(nbrecv,p)%omz*((j-0.5)*dy-anb(nbrecv,p)%y + boundfrontmyid)
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
                            if(sum1/sum2.gt.gamu(i,j,k).and.itype.eq.1) then
                               gamu(i,j,k) = sum1/sum2
                               upart(i,j,k) = anb(nbrecv,p)%u + &
                                    anb(nbrecv,p)%omy*((k-0.5)*dz-anb(nbrecv,p)%z) - &
                                    anb(nbrecv,p)%omz*((j-0.5)*dy-anb(nbrecv,p)%y + boundfrontmyid)
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
                         if (dist2 .lt. radin2.and.itype.eq.1) then
                            gamv(i,j,k) = 1.
                            vpart(i,j,k) = anb(nbrecv,p)%v + &
                                 anb(nbrecv,p)%omz*((i-0.5)*dx-anb(nbrecv,p)%x + boundleftmyid) - &
                                 anb(nbrecv,p)%omx*((k-0.5)*dz-anb(nbrecv,p)%z)
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
                            if(sum1/sum2.gt.gamv(i,j,k).and.itype.eq.1) then
                               gamv(i,j,k) = sum1/sum2
                               vpart(i,j,k) = anb(nbrecv,p)%v + &
                                    anb(nbrecv,p)%omz*((i-0.5)*dx-anb(nbrecv,p)%x + boundleftmyid) - &
                                    anb(nbrecv,p)%omx*((k-0.5)*dz-anb(nbrecv,p)%z)
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
                         if (dist2 .lt. radin2.and.itype.eq.1) then
                            gamw(i,j,k) = 1. 
                            wpart(i,j,k) = anb(nbrecv,p)%w + &
                                 anb(nbrecv,p)%omx*((j-0.5)*dy-anb(nbrecv,p)%y + boundfrontmyid) - &
                                 anb(nbrecv,p)%omy*((i-0.5)*dx-anb(nbrecv,p)%x + boundleftmyid)
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
                            if(sum1/sum2.gt.gamw(i,j,k).and.itype.eq.1) then
                               gamw(i,j,k) = sum1/sum2
                               wpart(i,j,k) = anb(nbrecv,p)%w + &
                                    anb(nbrecv,p)%omx*((j-0.5)*dy-anb(nbrecv,p)%y + boundfrontmyid) - &
                                    anb(nbrecv,p)%omy*((i-0.5)*dx-anb(nbrecv,p)%x + boundleftmyid)
                            endif
                         endif
                      enddo ! do i=
                   enddo ! do j=
                enddo ! do k=
                !
                ! pressure
                !
                do k=klow,khigh
                   coorzmin  = (k-1)*dz
                   coorzplus = (k)*dz
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
                         if (dist2 .lt. radin2.and.itype.eq.1) then
                            gamp(i,j,k) = 1. 
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
                            if(itype .eq. 1) gamp(i,j,k) = max(sum1/sum2,gamp(i,j,k)) 
                         endif
                         if (dist2 .lt. radin2.and.itype.eq.0) then
                            gamp(i,j,k) = 1. 
                            upart(i,j,k) = anb(nbrecv,p)%omx
                            vpart(i,j,k) = anb(nbrecv,p)%omy
                            wpart(i,j,k) = anb(nbrecv,p)%omz
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
                            if(sum1/sum2.gt.gamp(i,j,k).and.itype.eq.0) then
                               gamp(i,j,k) = sum1/sum2
                               upart(i,j,k) = anb(nbrecv,p)%omx
                               vpart(i,j,k) = anb(nbrecv,p)%omy
                               wpart(i,j,k) = anb(nbrecv,p)%omz
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
    gamp(:,:,:) = 1. - gamp(:,:,:)
    gamu(:,:,:) = 1. - gamu(:,:,:)
    gamv(:,:,:) = 1. - gamv(:,:,:)
    gamw(:,:,:) = 1. - gamw(:,:,:)

    call updthalos(gamp,1)
    call updthalos(gamu,1)
    call updthalos(gamv,1)
    call updthalos(gamw,1)
    call updthalos(gamp,2)
    call updthalos(gamu,2)
    call updthalos(gamv,2)
    call updthalos(gamw,2)
    call updthalos(upart,1)
    call updthalos(vpart,1)
    call updthalos(wpart,1)
    call updthalos(upart,2)
    call updthalos(vpart,2)
    call updthalos(wpart,2)
    !
    do j=0,j1
       do i=0,i1
          gamu(i,j,0)  = gamu(i,j,1)
          gamv(i,j,0)  = gamv(i,j,1)
          gamw(i,j,0)  = gamw(i,j,1)
          gamp(i,j,0)  = gamp(i,j,1)
          gamu(i,j,k1) = gamu(i,j,kmax)
          gamv(i,j,k1) = gamv(i,j,kmax)
          gamw(i,j,k1) = gamw(i,j,kmax)
          gamp(i,j,k1) = gamp(i,j,kmax)
       enddo
    enddo
    do j=0,j1
       do i=0,i1
          upart(i,j,0)    = -upart(i,j,1)     ! no-slip
          vpart(i,j,0)    = -vpart(i,j,1)     ! no-slip
          wpart(i,j,0)    = 0.                ! no-penetration
          upart(i,j,k1)   = -upart(i,j,kmax)  ! no-slip
          vpart(i,j,k1)   = -vpart(i,j,kmax)  ! no-slip
          wpart(i,j,kmax) = 0.                ! no-penetration
          wpart(i,j,k1)   = wpart(i,j,kmax-1) ! dw/dz=0 at wall (not used) 
       enddo
    enddo
    !
    !aux = sum(1.-gamp(1:imax,1:jmax,1:kmax))
    !call mpi_allreduce(aux,aux_all,1,mpi_real8,mpi_sum,MPI_COMM_WORLD,error)
    !if(myid.eq.0) print*,'Check Phase Indicator p:',np*volp,aux_all*dveul, &
    !                                               (np*volp-aux_all*dveul)/(np*volp)
    !aux = sum(1.-gamu(1:imax,1:jmax,1:kmax))
    !call mpi_allreduce(aux,aux_all,1,mpi_real8,mpi_sum,MPI_COMM_WORLD,error)
    !if(myid.eq.0) print*,'Check Phase Indicator u:',np*volp,aux_all*dveul, &
    !                                               (np*volp-aux_all*dveul)/(np*volp)
    !aux = sum(1.-gamv(1:imax,1:jmax,1:kmax))
    !call mpi_allreduce(aux,aux_all,1,mpi_real8,mpi_sum,MPI_COMM_WORLD,error)
    !if(myid.eq.0) print*,'Check Phase Indicator v:',np*volp,aux_all*dveul, &
    !                                               (np*volp-aux_all*dveul)/(np*volp)
    !aux = sum(1.-gamw(1:imax,1:jmax,1:kmax))
    !call mpi_allreduce(aux,aux_all,1,mpi_real8,mpi_sum,MPI_COMM_WORLD,error)
    !if(myid.eq.0) print*,'Check Phase Indicator w:',np*volp,aux_all*dveul, &
    !                                               (np*volp-aux_all*dveul)/(np*volp)
    !
    return
  end subroutine phase_indicator
  !
end module mod_phase_indicator
