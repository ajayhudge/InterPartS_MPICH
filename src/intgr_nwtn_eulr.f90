module mod_intgr_nwtn_eulr
  use mpi
  use mod_common
  use mod_common_mpi
  use mod_intgr_over_sphere
  use mod_collisions
  implicit none
  private
  public intgr_nwtn_eulr
contains
  !
  subroutine intgr_nwtn_eulr
    implicit none
    integer :: p,q,botw,topw,nb,nbsend,nbrec,iter,itermax
    real :: boundleftnb,boundrightnb,boundfrontnb,boundbacknb
    real, dimension(npmax) :: posx,posy,posz,postheta,posphi, &
         velx,vely,velz, &
         omgx,omgy,omgz
    integer, dimension(npmax) :: colrank
    integer :: sumcolrank,sumcolrank_all,rankmax
    real, dimension(2,1) ::maxerr,maxerr_all
    real :: err,maxerror
    !
    integer, dimension(0:8) :: pmax_nb
    integer, dimension(1:npmax,0:8) :: mslv_nb,newmaster_nb
    ! newmaster:
    !   > 0 : nr of neighbor that has become the new master of particle p
    !     0 : a) myid didn't contain particle number p, or
    !         b) myid remains master of this particle, or
    !         c) myid remains slave of this particle
    real, dimension(npmax) :: colflgx,colflgy,colflgz
    real :: deltax,deltay,deltaz,deltan,dist,nx,ny,nz
    !
    integer,parameter :: nprocs = dims(1)*dims(2)
    integer :: k,i
    type neighbour
       real :: x,y,z,theta,phi,u,v,w,omx,omy,omz
    end type neighbour
    type(neighbour), dimension(1:npmax,0:8) :: anb
    !type neighbour2
    !  real :: x,y,z
    !end type neighbour2
    !type(neighbour2), dimension(1:npmax,0:8) :: anb2
    integer :: tag
    integer :: nrrequests
    !integer :: arrayrequests(1:30)
    !integer :: arraystatuses(MPI_STATUS_SIZE,1:30)
!     integer :: arrayrequests(1:2)                      !original
!     integer :: arraystatuses(MPI_STATUS_SIZE,1:2)      !original
    integer, parameter :: MAX_REQ = 64                   !added
    integer :: arrayrequests(1:MAX_REQ)                  !added
    integer :: arraystatuses(MPI_STATUS_SIZE,1:MAX_REQ)  !added
    real :: ax,ay
    integer :: l
    real :: leftbound,rightbound,frontbound,backbound
    integer :: nbrec2
    integer, dimension(ndims) :: proccoords
    integer :: procrank
    integer :: counter !delete after checking things
    real :: eps
    integer :: count_mstr,count_slve,count_mstr_all,count_slve_loc
    integer :: idp,idp_nb,idq
    logical :: found_mstr
    integer :: qlast,qq
    character(len=3) :: rankpr
    real, dimension(3) :: arrdeltax,arrdeltay
    integer :: ideltax,ideltay
!    real :: isperiodx,isperiody
    real :: xfploc,yfploc,zfploc
    real :: isperiodx,isperiody
    real :: coorxfp,cooryfp,coorzfp
    logical :: isout
    integer :: ll
    !
    ! Initialization: new --> old
    !
    !$omp parallel default(none) &
    !$omp&shared(ap,op,pmax) &
    !$omp&private(p)
    !$omp do 
    do p=1,pmax
       if (ap(p)%mslv.gt.0) then
          op(p)%x = ap(p)%x
          op(p)%y = ap(p)%y
          op(p)%z = ap(p)%z
          op(p)%theta = ap(p)%theta
          op(p)%phi = ap(p)%phi
          op(p)%u = ap(p)%u
          op(p)%v = ap(p)%v
          op(p)%w = ap(p)%w
          op(p)%omx = ap(p)%omx
          op(p)%omy = ap(p)%omy
          op(p)%omz = ap(p)%omz
          op(p)%omtheta = ap(p)%omtheta
          op(p)%intu = ap(p)%intu
          op(p)%intv = ap(p)%intv
          op(p)%intw = ap(p)%intw
          op(p)%intomx = ap(p)%intomx
          op(p)%intomy = ap(p)%intomy
          op(p)%intomz = ap(p)%intomz
          op(p)%colfx = ap(p)%colfx
          op(p)%colfy = ap(p)%colfy
          op(p)%colfz = ap(p)%colfz
          op(p)%coltx = ap(p)%coltx
          op(p)%colty = ap(p)%colty
          op(p)%coltz = ap(p)%coltz
          op(p)%dx(:) = ap(p)%dx(:) 
          op(p)%dy(:) = ap(p)%dy(:)
          op(p)%dz(:) = ap(p)%dz(:)
          op(p)%dxt(:) = ap(p)%dxt(:)
          op(p)%dyt(:) = ap(p)%dyt(:)
          op(p)%dzt(:) = ap(p)%dzt(:)
          op(p)%dut(:) = ap(p)%dut(:)
          op(p)%dvt(:) = ap(p)%dvt(:)
          op(p)%dwt(:) = ap(p)%dwt(:)
       endif
    enddo
    !$omp end parallel
    !
    ! compute integral of linear and angular momentum over sphere
    !
    !call intgr_over_sphere(ap(:)%intu,ap(:)%intv,ap(:)%intw,2)
    !call intgr_over_sphere(ap(:)%intomx,ap(:)%intomy,ap(:)%intomz,3)
    call intgr_over_sphere(2)
    call intgr_over_sphere(3)
    !
    ! exchange data with neighbors: pmasterslave
    !
    pmax_nb(:) = 0
    pmax_nb(0)=pmax
    !$omp workshare
    mslv_nb(:,:) = 0
    mslv_nb(1:pmax,0) = ap(1:pmax)%mslv
    !$omp end workshare
    do nb=1,8
       nbsend = 4+nb
       nbrec  = nb
       if (nbsend .gt. 8) nbsend = nbsend - 8
       call MPI_SENDRECV(pmax_nb(0),1,MPI_INTEGER,neighbor(nbsend),1, &
            pmax_nb(nbrec),1,MPI_INTEGER,neighbor(nbrec),1, &
            comm_cart,status,error)
       call MPI_SENDRECV(mslv_nb(1,0),pmax_nb(0),MPI_INTEGER,neighbor(nbsend),2, &
            mslv_nb(1,nbrec),pmax_nb(nbrec),MPI_INTEGER,neighbor(nbrec),2, &
            comm_cart,status,error)
    enddo
    !
    ! begin iterative loop
    !
    iter = -1
    itermax = 1
    maxerror = 1
    do while (iter.lt.itermax)!.and.(sumcolrank_all+pmax*Nproc).ne.0.and.maxerror*dxi.gt.1.e-5
       iter = iter + 1
       !$omp workshare
       colrank(1:pmax) = -1 ! -1 means that particle p is not involved in a collision
       !$omp end workshare

       !$omp parallel default(none) &
       !$omp&shared(ap,anb,pmax) &
       !$omp&shared(posx,posy,posz,postheta,posphi,velx,vely,velz,omgx,omgy,omgz) &
       !$omp&private(p) &
       !$omp&firstprivate(iter)  
       !$omp do 
       do p=1,pmax
          if (ap(p)%mslv .gt. 0) then
             !myid is master of particle ap(p)%mslv
             posx(p) = ap(p)%x
             posy(p) = ap(p)%y
             posz(p) = ap(p)%z
             postheta(p) = ap(p)%theta
             posphi(p) = ap(p)%phi
             velx(p) = ap(p)%u
             vely(p) = ap(p)%v
             velz(p) = ap(p)%w
             omgx(p) = ap(p)%omx
             omgy(p) = ap(p)%omy
             omgz(p) = ap(p)%omz
          else
             posx(p)    = 0.
             posy(p)    = 0.
             posz(p)    = 0.
             postheta(p) = 0. 
             posphi(p)   = 0. 
             velx(p)    = 0.
             vely(p)    = 0.
             velz(p)    = 0.
             omgx(p)    = 0.
             omgy(p)    = 0.
             omgz(p)    = 0. ! I think I can simply put anb(p,0)%x here, check later
          endif
       enddo
       !
       ! exchange data with neighbors: posx,posy,posz
       !
       !$omp do  
       do p=1,pmax
          anb(p,0)%x = posx(p)
          anb(p,0)%y = posy(p)
          anb(p,0)%z = posz(p)
          anb(p,0)%theta = postheta(p) 
          anb(p,0)%phi   = posphi(p)
          anb(p,0)%u = velx(p)
          anb(p,0)%v = vely(p)
          anb(p,0)%w = velz(p)
          anb(p,0)%omx = omgx(p)
          anb(p,0)%omy = omgy(p)
          anb(p,0)%omz = omgz(p)
       enddo
       !$omp end parallel
       do nb=1,8
          nbsend = 4+nb
          nbrec  = nb
          if (nbsend .gt. 8) nbsend = nbsend - 8
          call MPI_SENDRECV(anb(1,0)%x,pmax_nb(0)*11,MPI_REAL8,neighbor(nbsend),3, &
               anb(1,nbrec)%x,pmax_nb(nbrec)*11,MPI_REAL8,neighbor(nbrec),3, &
               comm_cart,status,error)
          ! send x,y,z,u,v,w,omx,omy,omz -> 9*pmax contiguous info
          ! (see definition of type neighbor in the begining of the subroutine)
       enddo
       !
       ! recompute particle positions because of periodic b.c.'s in the x and y-direction
       !
       !$omp parallel default(none) &
       !$omp&shared(ap,anb,pmax,pmax_nb,mslv_nb,coords) &
       !$omp&private(p,nb,boundleftnb,boundbacknb,boundrightnb,boundfrontnb)
       nb=1
       !$omp do
       do p=1,pmax_nb(nb)
          if (mslv_nb(p,nb) .gt. 0) then
             boundleftnb  = (coords(1)+1)*lx/(1.*dims(1)) ! left  boundary of neighbor nb
             if (anb(p,nb)%x .lt. boundleftnb) then
                anb(p,nb)%x = anb(p,nb)%x + lx
             endif
          endif
       enddo
       nb=2
       !$omp do
       do p=1,pmax_nb(nb)
          if (mslv_nb(p,nb) .gt. 0) then
             boundleftnb  = (coords(1)+1)*lx/(1.*dims(1)) ! left  boundary of neighbor nb
             boundbacknb  = (coords(2))*ly/(1.*dims(2)) ! back  boundary of neighbor nb
             if (anb(p,nb)%x .lt. boundleftnb) then
                anb(p,nb)%x = anb(p,nb)%x + lx
             endif
             if (anb(p,nb)%y .gt. boundbacknb) then
                anb(p,nb)%y = anb(p,nb)%y - ly
             endif
          endif
       enddo
       nb=3
       !$omp do
       do p=1,pmax_nb(nb)
          if (mslv_nb(p,nb) .gt. 0) then
             boundbacknb  = (coords(2))*ly/(1.*dims(2)) ! back  boundary of neighbor nb
             if (anb(p,nb)%y .gt. boundbacknb) then
                anb(p,nb)%y = anb(p,nb)%y - ly
             endif
          endif
       enddo
       nb=4
       !$omp do
       do p=1,pmax_nb(nb)
          if (mslv_nb(p,nb) .gt. 0) then
             boundrightnb = (coords(1))*lx/(1.*dims(1)) ! right boundary of neighbor nb
             boundbacknb  = (coords(2))*ly/(1.*dims(2)) ! back  boundary of neighbor nb
             if (anb(p,nb)%x .gt. boundrightnb) then
                anb(p,nb)%x = anb(p,nb)%x - lx
             endif
             if (anb(p,nb)%y .gt. boundbacknb) then
                anb(p,nb)%y = anb(p,nb)%y - ly
             endif
          endif
       enddo
       nb=5
       !$omp do
       do p=1,pmax_nb(nb)
          if (mslv_nb(p,nb) .gt. 0) then
             boundrightnb = (coords(1))*lx/(1.*dims(1)) ! right boundary of neighbor nb
             if (anb(p,nb)%x .gt. boundrightnb) then
                anb(p,nb)%x = anb(p,nb)%x - lx
             endif
          endif
       enddo
       nb=6
       !$omp do
       do p=1,pmax_nb(nb)
          if (mslv_nb(p,nb) .gt. 0) then
             boundrightnb = (coords(1))*lx/(1.*dims(1)) ! right boundary of neighbor nb
             boundfrontnb = (coords(2)+1)*ly/(1.*dims(2)) ! front boundary of neighbor nb
             if (anb(p,nb)%x .gt. boundrightnb) then
                anb(p,nb)%x = anb(p,nb)%x - lx
             endif
             if (anb(p,nb)%y .lt. boundfrontnb) then
                anb(p,nb)%y = anb(p,nb)%y + ly
             endif
          endif
       enddo
       nb=7
       !$omp do
       do p=1,pmax_nb(nb)
          if (mslv_nb(p,nb) .gt. 0) then
             boundfrontnb = (coords(2)+1)*ly/(1.*dims(2)) ! front boundary of neighbor nb
             if (anb(p,nb)%y .lt. boundfrontnb) then
                anb(p,nb)%y = anb(p,nb)%y + ly
             endif
          endif
       enddo
       nb=8
       !$omp do
       do p=1,pmax_nb(nb)
          if (mslv_nb(p,nb) .gt. 0) then
             boundleftnb  = (coords(1)+1)*lx/(1.*dims(1)) ! left  boundary of neighbor nb
             boundfrontnb = (coords(2)+1)*ly/(1.*dims(2)) ! front boundary of neighbor nb
             if (anb(p,nb)%x .lt. boundleftnb) then
                anb(p,nb)%x = anb(p,nb)%x + lx
             endif
             if (anb(p,nb)%y .lt. boundfrontnb) then
                anb(p,nb)%y = anb(p,nb)%y + ly
             endif
          endif
       enddo
       !$omp end parallel
       !
       ! collision model
       !
       botw = np + 1 ! id of the bottom wall
       topw = np + 2 ! id of the top wall
       !
       !$omp parallel default(none) &
       !$omp&shared(ap,op,anb,pmax) &
       !$omp&shared(pmax_nb,mslv_nb) &
       !$omp&shared(colflgx,colflgy,colflgz,colrank) &
       !$omp&shared(itermax,myid,dt,rkiter,topw,botw) &
       !$omp&private(p,q,nb,dist,deltax,deltay,deltaz,deltan,eps) &
       !$omp&private(nx,ny,nz) &
       !$omp&private(idp,idq,qq,qlast) &
       !$omp&private(arrdeltax,arrdeltay) &
       !$omp&private(ideltax,ideltay) &
       !$omp&firstprivate(iter)
       !$omp do
       do p=1,pmax
          qlast = ap(p)%qmax
          colflgx(p) = 0.
          colflgy(p) = 0.
          colflgz(p) = 0.
          if(ap(p)%mslv.gt.0) then
             idp = ap(p)%mslv
             qlast = ap(p)%qmax
             ap(p)%colfx = 0.
             ap(p)%colfy = 0.
             ap(p)%colfz = 0.
             ap(p)%coltx = 0.
             ap(p)%colty = 0.
             ap(p)%coltz = 0.
             !
             ! (i) collision with other particles
             !
             do nb=0,8
                do q=1,pmax_nb(nb)
                   if((nb.ge.0).and.(idp.eq.mslv_nb(q,nb))) then
                      dist = 0. ! dummy operation
                   elseif(mslv_nb(q,nb).gt.0) then
                      idq = mslv_nb(q,nb)
                      !          deltax = anb(q,nb)%x - ap(p)%x
                      !          deltay = anb(q,nb)%y - ap(p)%y
                      deltaz = anb(q,nb)%z - ap(p)%z
                      arrdeltax(1) = anb(q,nb)%x - ap(p)%x
                      arrdeltax(2) = anb(q,nb)%x - lx - ap(p)%x
                      arrdeltax(3) = anb(q,nb)%x + lx - ap(p)%x
                      ideltax = minloc(abs(arrdeltax(1:3)),1)
                      deltax = arrdeltax(ideltax)
                      !
                      arrdeltay(1) = anb(q,nb)%y - ap(p)%y
                      arrdeltay(2) = anb(q,nb)%y - ly - ap(p)%y
                      arrdeltay(3) = anb(q,nb)%y + ly - ap(p)%y
                      ideltay = minloc(abs(arrdeltay(1:3)),1)
                      deltay = arrdeltay(ideltay)
                      !
                      dist = sqrt(deltax**2.+deltay**2.+deltaz**2.)
                      deltan = 2*radius-dist
                      nx = deltax/dist
                      ny = deltay/dist
                      nz = deltaz/dist ! computed twice (here and in the subroutine collisions)
                      eps = (dist-2.*radius)/radius
                      if((eps.lt.eps_ini_pp).and.(eps.gt.eps_cut_pp)) then
                         colrank(p) = myid
                         call lubrication(p,q,idp,idq,nx,ny,nz,eps, &
                              anb(q,nb)%u,anb(q,nb)%v,anb(q,nb)%w, &
                              anb(q,nb)%omx,anb(q,nb)%omy,anb(q,nb)%omz)
                      endif
                      qq = 0
                      !
                      ! Matching global id 'idq' with local id 'qq' of the particle in contact.
                      ! This has to be done because of the oblique collison model, which has 
                      ! memory: psi must be fixed at first conact and the tangential displacement
                      ! has to be integrated from the imminence of contact
                      !
                      do i=1,ap(p)%qmax
                         if(idq.eq.ap(p)%firstc(i)) then
                            !
                            ! ap(p)%firstc(i).eq.idq -> particles with ids 'idp' and 'idq' were in
                            ! contact in the previous substep
                            !
                            if(dist.lt.(2.*radius)) then 
                               !
                               ! particles are still in contact
                               !
                               qq = i
                               call collisions(p,q,qq,idp,idq,dist,deltax,deltay,deltaz, &
                                    anb(q,nb)%u,anb(q,nb)%v,anb(q,nb)%w, &
                                    anb(q,nb)%omx,anb(q,nb)%omy,anb(q,nb)%omz,dt)
                               !                if (abs(deltan*nx).gt.(colthr_pp)) colflgx(p) = 0. ! IBM force off
                               !                if (abs(deltan*ny).gt.(colthr_pp)) colflgy(p) = 0. ! IBM force off
                               !                if (abs(deltan*nz).gt.(colthr_pp)) colflgz(p) = 0. ! IBM force off
                               colrank(p) = myid
                            elseif(iter.eq.itermax) then
                               !
                               ! particles ceased to be in contact, remove local id and put the last
                               ! element in position 'i' to avoid 'gaps' in the local contact arrays
                               !
                               if(qlast.gt.1) then
                                  ap(p)%firstc(i) = ap(p)%firstc(qlast)
                                  ap(p)%dxt(i) = ap(p)%dxt(qlast)
                                  ap(p)%dyt(i) = ap(p)%dyt(qlast)
                                  ap(p)%dzt(i) = ap(p)%dzt(qlast)
                                  ap(p)%dut(i) = ap(p)%dut(qlast)
                                  ap(p)%dvt(i) = ap(p)%dvt(qlast)
                                  ap(p)%dwt(i) = ap(p)%dwt(qlast)
                                  op(p)%firstc(i) = op(p)%firstc(qlast)
                                  op(p)%dxt(i) = op(p)%dxt(qlast)
                                  op(p)%dyt(i) = op(p)%dyt(qlast)
                                  op(p)%dzt(i) = op(p)%dzt(qlast)
                                  op(p)%dut(i) = op(p)%dut(qlast)
                                  op(p)%dvt(i) = op(p)%dvt(qlast)
                                  op(p)%dwt(i) = op(p)%dwt(qlast)
                               endif
                               ap(p)%firstc(qlast) = 0
                               ap(p)%dx(qlast) = 0.
                               ap(p)%dy(qlast) = 0.
                               ap(p)%dz(qlast) = 0.
                               ap(p)%dxt(qlast) = 0.
                               ap(p)%dyt(qlast) = 0.
                               ap(p)%dzt(qlast) = 0.
                               ap(p)%dut(qlast) = 0.
                               ap(p)%dvt(qlast) = 0.
                               ap(p)%dwt(qlast) = 0.
                               op(p)%firstc(qlast) = 0
                               op(p)%dx(qlast) = 0.
                               op(p)%dy(qlast) = 0.
                               op(p)%dz(qlast) = 0.
                               op(p)%dxt(qlast) = 0.
                               op(p)%dyt(qlast) = 0.
                               op(p)%dzt(qlast) = 0.
                               op(p)%dut(qlast) = 0.
                               op(p)%dvt(qlast) = 0.
                               op(p)%dwt(qlast) = 0.
                               qlast = qlast - 1
                            endif
                         endif
                      enddo
                      if(qq.eq.0.and.dist.lt.(2.*radius)) then 
                         !
                         ! particle 'idp' in contact with particle 'idq' for the first time:
                         ! increase extent of local contact array and initialize it
                         !
                         qlast = qlast + 1
                         qq = qlast
                         ap(p)%dx(qlast) = 0.
                         ap(p)%dy(qlast) = 0.
                         ap(p)%dz(qlast) = 0.
                         ap(p)%dxt(qlast) = 0.
                         ap(p)%dyt(qlast) = 0.
                         ap(p)%dzt(qlast) = 0.
                         ap(p)%dut(qlast) = 0.
                         ap(p)%dvt(qlast) = 0.
                         ap(p)%dwt(qlast) = 0.
                         op(p)%dx(qlast) = 0.
                         op(p)%dy(qlast) = 0.
                         op(p)%dz(qlast) = 0.
                         op(p)%dxt(qlast) = 0.
                         op(p)%dyt(qlast) = 0.
                         op(p)%dzt(qlast) = 0.
                         op(p)%dut(qlast) = 0.
                         op(p)%dvt(qlast) = 0.
                         op(p)%dwt(qlast) = 0.
                         call collisions(p,q,qq,idp,idq,dist,deltax,deltay,deltaz, &
                              anb(q,nb)%u,anb(q,nb)%v,anb(q,nb)%w, &
                              anb(q,nb)%omx,anb(q,nb)%omy,anb(q,nb)%omz,dt)
                         !            if (abs(deltan*nx).gt.(colthr_pp)) colflgx(p) = 0. ! IBM force off
                         !            if (abs(deltan*ny).gt.(colthr_pp)) colflgy(p) = 0. ! IBM force off
                         !            if (abs(deltan*nz).gt.(colthr_pp)) colflgz(p) = 0. ! IBM force off
                         ap(p)%firstc(qlast) = idq
                         op(p)%firstc(qlast) = idq
                         colrank(p) = myid
                      endif
                      ap(p)%qmax = qlast ! size of local contact arrays updated
                   endif
                enddo
             enddo
             !
             ! (ii) collision with walls
             ! change velocities of the wall if the wall has a prescribed velocity (e.g. Couette flow)
             !
             do q=botw,topw
                idq = q
                deltax = 0.
                deltay = 0.
                deltaz = (q-botw)*lz - ap(p)%z
                dist = sqrt(deltax**2.+deltay**2.+deltaz**2.)
                deltan = radius-dist
                nx = deltax/dist
                ny = deltay/dist
                nz = deltaz/dist
                qq = 0
                eps = -deltan/radius 
                if((eps.lt.eps_ini_pw).and.(eps.gt.eps_cut_pw)) then
                   colrank(p) = myid
                   call lubrication(p,q,idp,idq,nx,ny,nz,eps, &
                        0.,0.,0.,0.,0.,0.)
                endif
                do i=1,ap(p)%qmax
                   if(idq.eq.ap(p)%firstc(i)) then
                      !
                      ! algorithm is analogous for particle-wall interactions, see comments
                      ! for particle-particle interactions 
                      !
                      qq = i
                      if(dist.lt.radius) then
                         call collisions(p,q,qq,idp,idq,dist,deltax,deltay,deltaz, &
                              0.,0.,0.,0.,0.,0.,dt)
                         if (dist.lt.(radius-colthr_pw)) colflgz(p) = 1. !IBM force off
                         colrank(p) = myid
                      elseif(iter.eq.itermax) then
                         if(qlast.gt.1) then
                            ap(p)%firstc(i) = ap(p)%firstc(qlast)
                            ap(p)%dxt(i) = ap(p)%dxt(qlast)
                            ap(p)%dyt(i) = ap(p)%dyt(qlast)
                            ap(p)%dzt(i) = ap(p)%dzt(qlast)
                            ap(p)%dut(i) = ap(p)%dut(qlast)
                            ap(p)%dvt(i) = ap(p)%dvt(qlast)
                            ap(p)%dwt(i) = ap(p)%dwt(qlast)
                            op(p)%firstc(i) = op(p)%firstc(qlast)
                            op(p)%dxt(i) = op(p)%dxt(qlast)
                            op(p)%dyt(i) = op(p)%dyt(qlast)
                            op(p)%dzt(i) = op(p)%dzt(qlast)
                            op(p)%dut(i) = op(p)%dut(qlast)
                            op(p)%dvt(i) = op(p)%dvt(qlast)
                            op(p)%dwt(i) = op(p)%dwt(qlast)
                         endif
                         ap(p)%firstc(qlast) = 0
                         ap(p)%dx(qlast) = 0.
                         ap(p)%dy(qlast) = 0.
                         ap(p)%dz(qlast) = 0.
                         ap(p)%dxt(qlast) = 0.
                         ap(p)%dyt(qlast) = 0.
                         ap(p)%dzt(qlast) = 0.
                         ap(p)%dut(qlast) = 0.
                         ap(p)%dvt(qlast) = 0.
                         ap(p)%dwt(qlast) = 0.
                         op(p)%firstc(qlast) = 0
                         op(p)%dx(qlast) = 0.
                         op(p)%dy(qlast) = 0.
                         op(p)%dz(qlast) = 0.
                         op(p)%dxt(qlast) = 0.
                         op(p)%dyt(qlast) = 0.
                         op(p)%dzt(qlast) = 0.
                         op(p)%dut(qlast) = 0.
                         op(p)%dvt(qlast) = 0.
                         op(p)%dwt(qlast) = 0.
                         qlast = qlast - 1
                      endif
                   endif
                enddo
                if((qq.eq.0).and.(dist.lt.radius)) then
                   qlast = qlast + 1
                   qq = qlast
                   ap(p)%dx(qlast) = 0.
                   ap(p)%dy(qlast) = 0.
                   ap(p)%dz(qlast) = 0.
                   ap(p)%dxt(qlast) = 0.
                   ap(p)%dyt(qlast) = 0.
                   ap(p)%dzt(qlast) = 0.
                   ap(p)%dut(qlast) = 0.
                   ap(p)%dvt(qlast) = 0.
                   ap(p)%dwt(qlast) = 0.
                   op(p)%dx(qlast) = 0.
                   op(p)%dy(qlast) = 0.
                   op(p)%dz(qlast) = 0.
                   op(p)%dxt(qlast) = 0.
                   op(p)%dyt(qlast) = 0.
                   op(p)%dzt(qlast) = 0.
                   op(p)%dut(qlast) = 0.
                   op(p)%dvt(qlast) = 0.
                   op(p)%dwt(qlast) = 0.
                   call collisions(p,q,qq,idp,idq,dist,deltax,deltay,deltaz, &
                        0.,0.,0.,0.,0.,0.,dt)
                   if ((dist).lt.(radius-colthr_pw)) colflgz(p) = 1. !switch off IBM force
                   ap(p)%firstc(qlast) = idq
                   op(p)%firstc(qlast) = idq
                   colrank(p) = myid
                endif
             enddo
             ap(p)%qmax = qlast
          endif
       enddo
       !
       !compute new particle positions
       !
!!$omp do schedule(dynamic)
       !$omp do 
       do p=1,pmax
          if (ap(p)%mslv.gt.0) then
             colflgx(p) = 0.
             colflgy(p) = 0.
             colflgz(p) = 0.
             ap(p)%u = op(p)%u + &  ! translational motion equation (SS)
                  (1-colflgx(p))*( &
                  (-1.)*dt*(ap(p)%fxltot)/(ap(p)%vol*ap(p)%ratiorho) + &
                  (ap(p)%intu-op(p)%intu)/(ap(p)%vol*ap(p)%ratiorho) &
                  ) + &
                  rkcoeffab(rkiter)*dt*gaccx*(1.-(1./ap(p)%ratiorho)) + &
                  rkcoeffab(rkiter)*0.5*dt*(ap(p)%colfx+op(p)%colfx)/(ap(p)%vol*ap(p)%ratiorho)
             ap(p)%x = op(p)%x + rkcoeffab(rkiter)*0.5*dt*(ap(p)%u+op(p)%u)
             ap(p)%v = op(p)%v + &
                  (1-colflgy(p))*( &
                  (-1.)*dt*(ap(p)%fyltot)/(ap(p)%vol*ap(p)%ratiorho) + &
                  (ap(p)%intv-op(p)%intv)/(ap(p)%vol*ap(p)%ratiorho) &
                  ) + &
                  rkcoeffab(rkiter)*dt*gaccy*(1.-(1./ap(p)%ratiorho)) + &
                  rkcoeffab(rkiter)*0.5*dt*(ap(p)%colfy+op(p)%colfy)/(ap(p)%vol*ap(p)%ratiorho)
             ap(p)%y = op(p)%y + rkcoeffab(rkiter)*0.5*dt*(ap(p)%v+op(p)%v)
             ap(p)%w = op(p)%w + &
                  (1-colflgz(p))*( &
                  (-1.)*dt*(ap(p)%fzltot)/(ap(p)%vol*ap(p)%ratiorho) + &
                  (ap(p)%intw-op(p)%intw)/(ap(p)%vol*ap(p)%ratiorho) &
                  ) + &
                  rkcoeffab(rkiter)*dt*gaccz*(1.-(1./ap(p)%ratiorho)) + &
                  rkcoeffab(rkiter)*0.5*dt*(ap(p)%colfz+op(p)%colfz)/(ap(p)%vol*ap(p)%ratiorho)
             ap(p)%z = op(p)%z + rkcoeffab(rkiter)*0.5*dt*(ap(p)%w+op(p)%w)
             !
             ! Rotational motion equations (SS) - with control switch
             if (particle_rotation) then
                ap(p)%omx = op(p)%omx + &
                     (-1.)*dt*ap(p)%torqxltot/(ap(p)%mominert*ap(p)%ratiorho) + &
                     (ap(p)%intomx-op(p)%intomx)/(ap(p)%mominert*ap(p)%ratiorho) + &
                     rkcoeffab(rkiter)*0.5*dt*(ap(p)%coltx+op(p)%coltx)/(ap(p)%mominert*ap(p)%ratiorho)
                ap(p)%omy = op(p)%omy + &
                     (-1.)*dt*ap(p)%torqyltot/(ap(p)%mominert*ap(p)%ratiorho) + &
                     (ap(p)%intomy-op(p)%intomy)/(ap(p)%mominert*ap(p)%ratiorho) + &
                     rkcoeffab(rkiter)*0.5*dt*(ap(p)%colty+op(p)%colty)/(ap(p)%mominert*ap(p)%ratiorho)
                ap(p)%omz = op(p)%omz + &
                     (-1.)*dt*ap(p)%torqzltot/(ap(p)%mominert*ap(p)%ratiorho) + &
                     (ap(p)%intomz-op(p)%intomz)/(ap(p)%mominert*ap(p)%ratiorho) + &
                     rkcoeffab(rkiter)*0.5*dt*(ap(p)%coltz+op(p)%coltz)/(ap(p)%mominert*ap(p)%ratiorho)
             else
                ! No rotation: set all angular velocities to zero
                ap(p)%omx = 0.0
                ap(p)%omy = 0.0
                ap(p)%omz = 0.0
             endif
             ap(p)%phi   = op(p)%phi + 0.5*rkcoeffab(rkiter)*dt*(ap(p)%omz+op(p)%omz)

             ap(p)%omtheta = (ap(p)%omy*cos(ap(p)%phi)) - & 
                  (ap(p)%omx*sin(ap(p)%phi))
             ap(p)%theta = op(p)%theta + rkcoeffab(rkiter)*0.5*dt*(ap(p)%omtheta+op(p)%omtheta)
          endif
       enddo
       !$omp end parallel
       !
       ! compute particle forces
       !
       call sumrk3
       !
       !Check whether a collision occured: (sumcolrank_all + pmax*Nproc) .ne. 0
       !
       sumcolrank = sum(colrank(1:npmax))
       call mpi_allreduce(sumcolrank,sumcolrank_all,1,mpi_integer,mpi_sum, &
            comm_cart,error)
       if ((sumcolrank_all+npmax*Nproc) .ne. 0) then !criterion for occurence of collisions 
          !THIS HAS TO BE CORRECTED, WILL DO LATER!!!!
          maxerror = 0.
          do p=1,pmax
             if (ap(p)%mslv .gt. 0) then
                err = sqrt( (posx(p)-ap(p)%x)**2. + (posy(p)-ap(p)%y)**2. + (posz(p)-ap(p)%z)**2. )
                if (err .gt. maxerror) maxerror = err
             endif
          enddo
          maxerr(1,1) = maxerror
          maxerr(2,1) = myid*1.
          call mpi_allreduce(maxerr,maxerr_all,1,mpi_2double_precision,mpi_maxloc, &
               comm_cart,error)
          maxerror = maxerr_all(1,1)
          rankmax  = nint(maxerr_all(2,1))
          if ( myid .eq. rankmax.and.iter.eq.1) then
             write(6,*) 'Iter., (max error)/dx = ',iter,maxerror*dxi
          endif
       endif
    enddo ! do while
    !
    ! correction for periodic b.c.'s
    !
    !$omp parallel default(none) &
    !$omp&shared(ap,pmax) &
    !$omp&private(p)
    !$omp do
    do p=1,pmax
       if (ap(p)%mslv .gt. 0) then
          if (ap(p)%x .gt. lx) ap(p)%x = ap(p)%x-lx
          if (ap(p)%x .lt. 0.) ap(p)%x = ap(p)%x+lx
          if (ap(p)%y .gt. ly) ap(p)%y = ap(p)%y-ly
          if (ap(p)%y .lt. 0.) ap(p)%y = ap(p)%y+ly
       endif
    enddo
    !$omp end parallel
    !
    ! particle positions were updated. Now check if there are new masters
    !
    do p=1,pmax
       newmaster_nb(p,0) = 0
       if (ap(p)%mslv.gt.0) then
          ! myid was master of particle ap(p)%mslv at previous time step n
          ax = 0.5
          ay = 0.5
          if (ap(p)%x.eq.lx) ax = 0.51
          if (ap(p)%x.eq.0.) ax = 0.49
          if (ap(p)%y.eq.ly) ay = 0.51
          if (ap(p)%y.eq.0.) ay = 0.49
          proccoords(1) = nint( dims(1)*ap(p)%x/lx - ax )
          proccoords(2) = nint( dims(2)*ap(p)%y/ly - ay ) 
          call MPI_CART_RANK(comm_cart,proccoords,procrank,error)
          if (procrank .ne. myid) then
             ! particle ap(p)%mslv has a new master at time step n+1
             do nb=1,8
                if (procrank .eq. neighbor(nb)) then
                   newmaster_nb(p,0) = nb
                endif
             enddo
          endif
       endif
    enddo
    !
    ! exchange data
    !
    do nb=1,8
       nbsend = nb
       nbrec  = nb+4
       if (nbrec .gt. 8) nbrec = nbrec-8
       call MPI_SENDRECV(newmaster_nb(1,0),pmax_nb(0),MPI_INTEGER,neighbor(nbsend),1, &
            newmaster_nb(1,nbrec),pmax_nb(nbrec),MPI_INTEGER,neighbor(nbrec),1, &
            comm_cart,status,error)
    enddo
    !
    do p=1,pmax
       nrrequests = 0
       if (newmaster_nb(p,0).gt.0) then
          nbsend = newmaster_nb(p,0)
          tag = ap(p)%mslv!*10+nbsend
          nrrequests = nrrequests + 1
          call MPI_ISEND(ap(p)%x,send_real,MPI_REAL8,neighbor(nbsend),tag,comm_cart,arrayrequests((nrrequests-1)*3+1),error)
          call MPI_ISEND(ap(p)%qmax,send_int,MPI_INTEGER,neighbor(nbsend),tag+np,comm_cart,arrayrequests((nrrequests-1)*3+2),error)
          ap(p)%mslv = -ap(p)%mslv ! master is a slave now
          !
          ! NEW
          !
          call MPI_ISEND(rkp(p)%fxltot,19,MPI_REAL8,neighbor(nbsend),tag+2*np,comm_cart,arrayrequests((nrrequests-1)*3+3),error)
       endif
       if(mslv_nb(p,0).lt.0) then
          idp = -ap(p)%mslv
          do nbrec = 1,8
             nbrec2 = nbrec + 4
             if(nbrec2 .gt. 8) nbrec2 = nbrec2-8
             do i=1,pmax_nb(nbrec)
                idp_nb = mslv_nb(i,nbrec)
                if(newmaster_nb(i,nbrec) .eq. nbrec2.and.idp.eq.idp_nb) then
                   ap(p)%mslv = -ap(p)%mslv ! slave became a master
                   nrrequests = nrrequests + 1
                   tag = ap(p)%mslv!*10+nbsend
                   call MPI_IRECV(ap(p)%x,send_real,MPI_REAL8,neighbor(nbrec),tag,comm_cart,arrayrequests((nrrequests-1)*3+1),error)
                   call MPI_IRECV(ap(p)%qmax,send_int,MPI_INTEGER,neighbor(nbrec),tag+np,comm_cart,arrayrequests((nrrequests-1)*3+2),error)
                   !
                   ! NEW
                   !
                   call MPI_IRECV(rkp(p)%fxltot,19,MPI_REAL8,neighbor(nbrec),tag+2*np,comm_cart,arrayrequests((nrrequests-1)*3+3),error)
                endif
             enddo
          enddo
       endif
       nrrequests = nrrequests*3
       call MPI_WAITALL(nrrequests,arrayrequests,arraystatuses,error)
    enddo
    !
    ! Masters are known now: ap(p)%mslv > 0.  
    ! Next step: determine slaves and master/slave neighbors
    !
    forall(p=1:pmax,k=0:8) mslv_nb(p,k) = 0
    mslv_nb(1:pmax,0) = ap(1:pmax)%mslv  ! new masters, but new slaves not all determined yet!
    ! process might remain master, but slave might lose particle
    !$omp workshare
    anb(1:pmax,1:8)%x     = 0.
    anb(1:pmax,1:8)%y     = 0.
    anb(1:pmax,1:8)%z     = 0.
    anb(1:pmax,1:8)%theta = 0.
    anb(1:pmax,1:8)%phi   = 0.
    anb(1:pmax,1:8)%u     = 0.
    anb(1:pmax,1:8)%v     = 0.
    anb(1:pmax,1:8)%w     = 0.
    anb(1:pmax,1:8)%omx   = 0.
    anb(1:pmax,1:8)%omy   = 0.
    anb(1:pmax,1:8)%omz   = 0.
    anb(1:pmax,0)%x = ap(1:pmax)%x
    anb(1:pmax,0)%y = ap(1:pmax)%y
    anb(1:pmax,0)%z = ap(1:pmax)%z
    anb(1:pmax,0)%theta = ap(1:pmax)%theta
    anb(1:pmax,0)%phi = ap(1:pmax)%phi
    anb(1:pmax,0)%u = ap(1:pmax)%u
    anb(1:pmax,0)%v = ap(1:pmax)%v
    anb(1:pmax,0)%w = ap(1:pmax)%w
    anb(1:pmax,0)%omx = ap(1:pmax)%omx
    anb(1:pmax,0)%omy = ap(1:pmax)%omy
    anb(1:pmax,0)%omz = ap(1:pmax)%omz
    !$omp end workshare
    do nb=1,8
       nbsend = 4+nb
       if (nbsend .gt. 8) nbsend = nbsend - 8
       nbrec  = nb
       !  nbsend = nb
       !  nbrec  = nb+4
       !  if (nbrec .gt. 8) nbrec = nbrec-8
       call MPI_SENDRECV(mslv_nb(1,0),pmax_nb(0),MPI_INTEGER,neighbor(nbsend),1, &
            mslv_nb(1,nbrec),pmax_nb(nbrec),MPI_INTEGER,neighbor(nbrec),1, &
            comm_cart,status,error)
       call MPI_SENDRECV(anb(1,0)%x,pmax_nb(0)*11,MPI_REAL8,neighbor(nbsend),2, &
            anb(1,nbrec)%x,pmax_nb(nbrec)*11,MPI_REAL8,neighbor(nbrec),2, &
            comm_cart,status,error)
       ! send x,y,z -> 3*pmax contiguous info
       ! (see definition of type neighbor2 in the begining of the subroutine)
    enddo
    !
    ! recompute particle positions because of periodic b.c.'s in the x and y-direction
    !
    !$omp parallel default(none) &
    !$omp&shared(ap,anb,coords,pmax_nb,mslv_nb)  &
    !$omp&private(p,nb,boundleftnb,boundbacknb,boundrightnb,boundfrontnb) 
    nb=1
    !$omp do
    do p=1,pmax_nb(nb)
       if (mslv_nb(p,nb) .gt. 0) then
          boundleftnb  = (coords(1)+1)*lx/(1.*dims(1)) ! left  boundary of neighbor nb
          if (anb(p,nb)%x .lt. boundleftnb) then
             anb(p,nb)%x = anb(p,nb)%x + lx
          endif
       endif
    enddo
    nb=2
    !$omp do
    do p=1,pmax_nb(nb)
       if (mslv_nb(p,nb) .gt. 0) then
          boundleftnb  = (coords(1)+1)*lx/(1.*dims(1)) ! left  boundary of neighbor nb
          boundbacknb  = (coords(2))*ly/(1.*dims(2)) ! back  boundary of neighbor nb
          if (anb(p,nb)%x .lt. boundleftnb) then
             anb(p,nb)%x = anb(p,nb)%x + lx
          endif
          if (anb(p,nb)%y .gt. boundbacknb) then
             anb(p,nb)%y = anb(p,nb)%y - ly
          endif
       endif
    enddo
    nb=3
    !$omp do
    do p=1,pmax_nb(nb)
       if (mslv_nb(p,nb) .gt. 0) then
          boundbacknb  = (coords(2))*ly/(1.*dims(2)) ! back  boundary of neighbor nb
          if (anb(p,nb)%y .gt. boundbacknb) then
             anb(p,nb)%y = anb(p,nb)%y - ly
          endif
       endif
    enddo
    nb=4
    !$omp do
    do p=1,pmax_nb(nb)
       if (mslv_nb(p,nb) .gt. 0) then
          boundrightnb = (coords(1))*lx/(1.*dims(1)) ! right boundary of neighbor nb
          boundbacknb  = (coords(2))*ly/(1.*dims(2)) ! back  boundary of neighbor nb
          if (anb(p,nb)%x .gt. boundrightnb) then
             anb(p,nb)%x = anb(p,nb)%x - lx
          endif
          if (anb(p,nb)%y .gt. boundbacknb) then
             anb(p,nb)%y = anb(p,nb)%y - ly
          endif
       endif
    enddo
    nb=5
    !$omp do
    do p=1,pmax_nb(nb)
       if (mslv_nb(p,nb) .gt. 0) then
          boundrightnb = (coords(1))*lx/(1.*dims(1)) ! right boundary of neighbor nb
          if (anb(p,nb)%x .gt. boundrightnb) then
             anb(p,nb)%x = anb(p,nb)%x - lx
          endif
       endif
    enddo
    nb=6
    !$omp do
    do p=1,pmax_nb(nb)
       if (mslv_nb(p,nb) .gt. 0) then
          boundrightnb = (coords(1))*lx/(1.*dims(1)) ! right boundary of neighbor nb
          boundfrontnb = (coords(2)+1)*ly/(1.*dims(2)) ! front boundary of neighbor nb
          if (anb(p,nb)%x .gt. boundrightnb) then
             anb(p,nb)%x = anb(p,nb)%x - lx
          endif
          if (anb(p,nb)%y .lt. boundfrontnb) then
             anb(p,nb)%y = anb(p,nb)%y + ly
          endif
       endif
    enddo
    nb=7
    !$omp do
    do p=1,pmax_nb(nb)
       if (mslv_nb(p,nb) .gt. 0) then
          boundfrontnb = (coords(2)+1)*ly/(1.*dims(2)) ! front boundary of neighbor nb
          if (anb(p,nb)%y .lt. boundfrontnb) then
             anb(p,nb)%y = anb(p,nb)%y + ly
          endif
       endif
    enddo
    nb=8
    !$omp do
    do p=1,pmax_nb(nb)
       if (mslv_nb(p,nb) .gt. 0) then
          boundleftnb  = (coords(1)+1)*lx/(1.*dims(1)) ! left  boundary of neighbor nb
          boundfrontnb = (coords(2)+1)*ly/(1.*dims(2)) ! front boundary of neighbor nb
          if (anb(p,nb)%x .lt. boundleftnb) then
             anb(p,nb)%x = anb(p,nb)%x + lx
          endif
          if (anb(p,nb)%y .lt. boundfrontnb) then
             anb(p,nb)%y = anb(p,nb)%y + ly
          endif
       endif
    enddo
    !$omp end parallel
    !
    ! important info that should be kept:
    !
    !$omp workshare
    sp(1:pmax)%mslv = ap(1:pmax)%mslv
    sp(1:pmax)%x = ap(1:pmax)%x
    sp(1:pmax)%y = ap(1:pmax)%y
    sp(1:pmax)%z = ap(1:pmax)%z
    sp(1:pmax)%theta = ap(1:pmax)%theta
    sp(1:pmax)%phi = ap(1:pmax)%phi
    sp(1:pmax)%u = ap(1:pmax)%u
    sp(1:pmax)%v = ap(1:pmax)%v
    sp(1:pmax)%w = ap(1:pmax)%w
    sp(1:pmax)%omx = ap(1:pmax)%omx
    sp(1:pmax)%omy = ap(1:pmax)%omy
    sp(1:pmax)%omz = ap(1:pmax)%omz
    sp(1:pmax)%omtheta = ap(1:pmax)%omtheta
    sp(1:pmax)%intu = ap(1:pmax)%intu
    sp(1:pmax)%intv = ap(1:pmax)%intv
    sp(1:pmax)%intw = ap(1:pmax)%intw
    sp(1:pmax)%intomx = ap(1:pmax)%intomx
    sp(1:pmax)%intomy = ap(1:pmax)%intomy
    sp(1:pmax)%intomz = ap(1:pmax)%intomz
    sp(1:pmax)%colfx = ap(1:pmax)%colfx
    sp(1:pmax)%colfy = ap(1:pmax)%colfy
    sp(1:pmax)%colfz = ap(1:pmax)%colfz
    sp(1:pmax)%coltx = ap(1:pmax)%coltx
    sp(1:pmax)%colty = ap(1:pmax)%colty
    sp(1:pmax)%coltz = ap(1:pmax)%coltz
    sp(1:pmax)%qmax = ap(1:pmax)%qmax
    forall(q=1:nqmax)
       sp(1:pmax)%dx(q) = ap(1:pmax)%dx(q)
       sp(1:pmax)%dy(q) = ap(1:pmax)%dy(q)
       sp(1:pmax)%dz(q) = ap(1:pmax)%dz(q)
       sp(1:pmax)%dxt(q) = ap(1:pmax)%dxt(q)
       sp(1:pmax)%dyt(q) = ap(1:pmax)%dyt(q)
       sp(1:pmax)%dzt(q) = ap(1:pmax)%dzt(q)
       sp(1:pmax)%firstc(q) = ap(1:pmax)%firstc(q)
    end forall
    !
    ! clear structure ap for re-ordering:
    !
    ap(1:npmax)%mslv = 0
    ap(1:npmax)%x = 0.
    ap(1:npmax)%y = 0.
    ap(1:npmax)%z = 0.
    ap(1:npmax)%theta = 0.
    ap(1:npmax)%phi = 0.
    ap(1:npmax)%u = 0.
    ap(1:npmax)%v = 0.
    ap(1:npmax)%w = 0.
    ap(1:npmax)%omx = 0.
    ap(1:npmax)%omy = 0.
    ap(1:npmax)%omz = 0.
    ap(1:npmax)%omtheta = 0.
    ap(1:npmax)%intu = 0.
    ap(1:npmax)%intv = 0.
    ap(1:npmax)%intw = 0.
    ap(1:npmax)%intomx = 0.
    ap(1:npmax)%intomy = 0.
    ap(1:npmax)%intomz = 0.
    ap(1:npmax)%colfx = 0.
    ap(1:npmax)%colfy = 0.
    ap(1:npmax)%colfz = 0.
    ap(1:npmax)%coltx = 0.
    ap(1:npmax)%colty = 0.
    ap(1:npmax)%coltz = 0.
    forall (p=1:npmax)
       ap(p)%dx(:) = 0.
       ap(p)%dy(:) = 0.
       ap(p)%dz(:) = 0.
       ap(p)%dxt(:) = 0.
       ap(p)%dyt(:) = 0.
       ap(p)%dzt(:) = 0.
       ap(p)%firstc(:) = 0
       ap(p)%xfp(:) = 0.
       ap(p)%yfp(:) = 0.
       ap(p)%zfp(:) = 0.
       ap(p)%nb(:) = 0
    end forall
    !$omp end workshare
    call transfer_sumrk3
    !
    leftbound   = (coords(1)  )*lx/(1.*dims(1)) ! left  boundary of process myid
    rightbound  = (coords(1)+1)*lx/(1.*dims(1)) ! right boundary of process myid
    frontbound  = (coords(2)  )*ly/(1.*dims(2)) ! front boundary of process myid
    backbound   = (coords(2)+1)*ly/(1.*dims(2)) ! back  boundary of process myid
    !
    i = 0
    count_mstr = 0
    count_slve = 0
    do idp = 1,np
       found_mstr = .false.
       call binsearch(idp,mslv_nb(1:pmax_nb(0),0),pmax_nb(0),found_mstr,p)
       !  do p=1,pmax
       !    if(sp(p)%mslv.eq.idp) then
       !      found_mstr = .true.
       if(found_mstr) then
          i = i + 1
          ap(i)%mslv = sp(p)%mslv
          ap(i)%x = sp(p)%x
          ap(i)%y = sp(p)%y
          ap(i)%z = sp(p)%z
          ap(i)%theta = sp(p)%theta
          ap(i)%phi = sp(p)%phi
          ap(i)%u = sp(p)%u
          ap(i)%v = sp(p)%v
          ap(i)%w = sp(p)%w
          ap(i)%omx = sp(p)%omx
          ap(i)%omy = sp(p)%omy
          ap(i)%omz = sp(p)%omz
          ap(i)%omtheta = sp(p)%omtheta
          ap(i)%intu = sp(p)%intu
          ap(i)%intv = sp(p)%intv
          ap(i)%intw = sp(p)%intw
          ap(i)%intomx = sp(p)%intomx
          ap(i)%intomy = sp(p)%intomy
          ap(i)%intomz = sp(p)%intomz
          ap(i)%colfx = sp(p)%colfx
          ap(i)%colfy = sp(p)%colfy
          ap(i)%colfz = sp(p)%colfz
          ap(i)%coltx = sp(p)%coltx
          ap(i)%colty = sp(p)%colty
          ap(i)%coltz = sp(p)%coltz
          ap(i)%qmax = sp(p)%qmax
          ap(i)%dx(1:nqmax) = sp(p)%dx(1:nqmax)
          ap(i)%dy(1:nqmax) = sp(p)%dy(1:nqmax)
          ap(i)%dz(1:nqmax) = sp(p)%dz(1:nqmax)
          ap(i)%dxt(1:nqmax) = sp(p)%dxt(1:nqmax)
          ap(i)%dyt(1:nqmax) = sp(p)%dyt(1:nqmax)
          ap(i)%dzt(1:nqmax) = sp(p)%dzt(1:nqmax)
          ap(i)%firstc(1:nqmax) = sp(p)%firstc(1:nqmax)
          call sumrk3_newmaster(i,p)
          count_mstr = count_mstr + 1
          ! neighbor 1
          if ( sp(p)%x .gt. (rightbound-(radius+offset)) ) then
             ap(i)%nb(1) = 1 ! neighbor 1 is slave of particle ap(i)%mslv
          endif
          ! neighbor 2
          dist = sqrt( (rightbound-sp(p)%x)**2. + (frontbound-sp(p)%y)**2. )
          if ( abs(dist) .lt. (radius+offset) ) then
             ap(i)%nb(2) = 1 ! neighbor 2 is slave of particle ap(i)%mslv
          endif
          ! neighbor 3
          if ( sp(p)%y .lt. (frontbound+(radius+offset)) ) then
             ap(i)%nb(3) = 1 ! neighbor 3 is slave of particle ap(i)%mslv
          endif
          ! neighbor 4
          dist = sqrt( (leftbound-sp(p)%x)**2. + (frontbound-sp(p)%y)**2. )
          if ( abs(dist) .lt. (radius+offset) ) then
             ap(i)%nb(4) = 1 ! neighbor 4 is slave of particle ap(i)%mslv
          endif
          ! neighbor 5
          if ( sp(p)%x .lt. (leftbound+(radius+offset)) ) then
             ap(i)%nb(5) = 1 ! neighbor 5 is slave of particle ap(i)%mslv
          endif
          ! neighbor 6
          dist = sqrt( (leftbound-sp(p)%x)**2. + (backbound-sp(p)%y)**2. )
          if ( abs(dist) .lt. (radius+offset) ) then
             ap(i)%nb(6) = 1 ! neighbor 6 is slave of particle ap(i)%mslv
          endif
          ! neighbor 7
          if ( sp(p)%y .gt. (backbound-(radius+offset)) ) then
             ap(i)%nb(7) = 1 ! neighbor 7 is slave of particle ap(i)%mslv
          endif
          ! neighbor 8
          dist = sqrt( (rightbound-sp(p)%x)**2. + (backbound-sp(p)%y)**2. )
          if ( abs(dist) .lt. (radius+offset) ) then
             ap(i)%nb(8) = 1 ! neighbor 8 is slave of particle ap(i)%mslv
          endif
          !      call sumrk3(i,p)
          !    endif
          !  enddo
          ! if(.NOT.found_mstr) then
       else
          count_slve_loc = 0
          nb=5
          !    do k=1,pmax_nb(nb)
          !      if ( (mslv_nb(k,nb).eq.idp) ) then
          ! neighbor nb is master of particle mslv_nb(p,nb)
          call binsearch(idp,mslv_nb(1:pmax_nb(nb),nb),pmax_nb(nb),found_mstr,k)
          if(found_mstr) then
             dist = sqrt( (leftbound-anb(k,nb)%x)**2. )
             if ( dist .lt. (radius+offset) ) then
                if(count_slve_loc.eq. 0 ) i = i+1
                ap(i)%mslv  = -mslv_nb(k,nb) ! myid is slave of particle mslv_nb(k,nb)
                ap(i)%nb(nb) = 1             ! neighbor nb of myid is particle's master
                count_slve_loc = count_slve_loc + 1
                isperiodx = 0.
                isperiody = 0.
                if (anb(k,nb)%x.lt.0.) isperiodx =  1.
                if (anb(k,nb)%x.gt.lx) isperiodx = -1.
                if (anb(k,nb)%y.lt.0.) isperiody =  1.
                if (anb(k,nb)%y.gt.ly) isperiody = -1.
                ap(i)%x = anb(k,nb)%x+isperiodx*lx
                ap(i)%y = anb(k,nb)%y+isperiody*ly
                ap(i)%z = anb(k,nb)%z
                ap(i)%theta = anb(k,nb)%theta
                ap(i)%phi   = anb(k,nb)%phi
                ap(i)%u     = anb(k,nb)%u
                ap(i)%v     = anb(k,nb)%v
                ap(i)%w     = anb(k,nb)%w
                ap(i)%omx   = anb(k,nb)%omx
                ap(i)%omy   = anb(k,nb)%omy
                ap(i)%omz   = anb(k,nb)%omz
             endif
          endif
          !    enddo
          nb=6
          !    do k=1,pmax_nb(nb)
          !      if ( (mslv_nb(k,nb).eq.idp) ) then
          ! neighbor nb is master of particle mslv_nb(p,nb)
          call binsearch(idp,mslv_nb(1:pmax_nb(nb),nb),pmax_nb(nb),found_mstr,k)
          if(found_mstr) then
             dist = sqrt( (leftbound-anb(k,nb)%x)**2. + (backbound-anb(k,nb)%y)**2. )
             if ( dist .lt. (radius+offset) ) then
                if(count_slve_loc.eq. 0 ) i = i+1
                ap(i)%mslv  = -mslv_nb(k,nb) ! myid is slave of particle mslv_nb(k,nb)
                ap(i)%nb(nb) = 1             ! neighbor nb of myid is particle's master
                count_slve_loc = count_slve_loc + 1
                isperiodx = 0.
                isperiody = 0.
                if (anb(k,nb)%x.lt.0.) isperiodx =  1.
                if (anb(k,nb)%x.gt.lx) isperiodx = -1.
                if (anb(k,nb)%y.lt.0.) isperiody =  1.
                if (anb(k,nb)%y.gt.ly) isperiody = -1.
                ap(i)%x = anb(k,nb)%x+isperiodx*lx
                ap(i)%y = anb(k,nb)%y+isperiody*ly
                ap(i)%z = anb(k,nb)%z
                ap(i)%theta = anb(k,nb)%theta
                ap(i)%phi   = anb(k,nb)%phi
                ap(i)%u     = anb(k,nb)%u
                ap(i)%v     = anb(k,nb)%v
                ap(i)%w     = anb(k,nb)%w
                ap(i)%omx   = anb(k,nb)%omx
                ap(i)%omy   = anb(k,nb)%omy
                ap(i)%omz   = anb(k,nb)%omz
             endif
          endif
          !    enddo
          nb=7
          !    do k=1,pmax_nb(nb)
          !      if ( (mslv_nb(k,nb).eq.idp) ) then
          ! neighbor nb is master of particle mslv_nb(p,nb)
          call binsearch(idp,mslv_nb(1:pmax_nb(nb),nb),pmax_nb(nb),found_mstr,k)
          if(found_mstr) then
             dist = sqrt( (backbound-anb(k,nb)%y)**2. )
             if ( dist .lt. (radius+offset) ) then
                if(count_slve_loc.eq. 0 ) i = i+1
                ap(i)%mslv  = -mslv_nb(k,nb) ! myid is slave of particle mslv_nb(k,nb)
                ap(i)%nb(nb) = 1             ! neighbor nb of myid is particle's master
                count_slve_loc = count_slve_loc + 1
                isperiodx = 0.
                isperiody = 0.
                if (anb(k,nb)%x.lt.0.) isperiodx =  1.
                if (anb(k,nb)%x.gt.lx) isperiodx = -1.
                if (anb(k,nb)%y.lt.0.) isperiody =  1.
                if (anb(k,nb)%y.gt.ly) isperiody = -1.
                ap(i)%x = anb(k,nb)%x+isperiodx*lx
                ap(i)%y = anb(k,nb)%y+isperiody*ly
                ap(i)%z = anb(k,nb)%z
                ap(i)%theta = anb(k,nb)%theta
                ap(i)%phi   = anb(k,nb)%phi
                ap(i)%u     = anb(k,nb)%u
                ap(i)%v     = anb(k,nb)%v
                ap(i)%w     = anb(k,nb)%w
                ap(i)%omx   = anb(k,nb)%omx
                ap(i)%omy   = anb(k,nb)%omy
                ap(i)%omz   = anb(k,nb)%omz
             endif
          endif
          !    enddo
          nb=8
          !    do k=1,pmax_nb(nb)
          !      if ( (mslv_nb(k,nb).eq.idp) ) then
          !neighbor nb is master of particle mslv_nb(p,nb)
          call binsearch(idp,mslv_nb(1:pmax_nb(nb),nb),pmax_nb(nb),found_mstr,k)
          if(found_mstr) then
             dist = sqrt( (rightbound-anb(k,nb)%x)**2. + (backbound-anb(k,nb)%y)**2. )
             if ( dist .lt. (radius+offset) ) then
                if(count_slve_loc.eq. 0 ) i = i+1
                ap(i)%mslv  = -mslv_nb(k,nb) ! myid is slave of particle mslv_nb(k,nb)
                ap(i)%nb(nb) = 1             ! neighbor nb of myid is particle's master
                count_slve_loc = count_slve_loc + 1
                isperiodx = 0.
                isperiody = 0.
                if (anb(k,nb)%x.lt.0.) isperiodx =  1.
                if (anb(k,nb)%x.gt.lx) isperiodx = -1.
                if (anb(k,nb)%y.lt.0.) isperiody =  1.
                if (anb(k,nb)%y.gt.ly) isperiody = -1.
                ap(i)%x = anb(k,nb)%x+isperiodx*lx
                ap(i)%y = anb(k,nb)%y+isperiody*ly
                ap(i)%z = anb(k,nb)%z
                ap(i)%theta = anb(k,nb)%theta
                ap(i)%phi   = anb(k,nb)%phi
                ap(i)%u     = anb(k,nb)%u
                ap(i)%v     = anb(k,nb)%v
                ap(i)%w     = anb(k,nb)%w
                ap(i)%omx   = anb(k,nb)%omx
                ap(i)%omy   = anb(k,nb)%omy
                ap(i)%omz   = anb(k,nb)%omz
             endif
          endif
          !    enddo
          nb=1
          !    do k=1,pmax_nb(nb)
          !      if ( (mslv_nb(k,nb).eq.idp) ) then
          ! neighbor nb is master of particle mslv_nb(p,nb)
          call binsearch(idp,mslv_nb(1:pmax_nb(nb),nb),pmax_nb(nb),found_mstr,k)
          if(found_mstr) then
             dist = sqrt( (rightbound-anb(k,nb)%x)**2. )
             if ( dist .lt. (radius+offset) ) then
                if(count_slve_loc.eq. 0 ) i = i+1
                ap(i)%mslv  = -mslv_nb(k,nb) ! myid is slave of particle mslv_nb(k,nb)
                ap(i)%nb(nb) = 1             ! neighbor nb of myid is particle's master
                count_slve_loc = count_slve_loc + 1
                isperiodx = 0.
                isperiody = 0.
                if (anb(k,nb)%x.lt.0.) isperiodx =  1.
                if (anb(k,nb)%x.gt.lx) isperiodx = -1.
                if (anb(k,nb)%y.lt.0.) isperiody =  1.
                if (anb(k,nb)%y.gt.ly) isperiody = -1.
                ap(i)%x = anb(k,nb)%x+isperiodx*lx
                ap(i)%y = anb(k,nb)%y+isperiody*ly
                ap(i)%z = anb(k,nb)%z
                ap(i)%theta = anb(k,nb)%theta
                ap(i)%phi   = anb(k,nb)%phi
                ap(i)%u     = anb(k,nb)%u
                ap(i)%v     = anb(k,nb)%v
                ap(i)%w     = anb(k,nb)%w
                ap(i)%omx   = anb(k,nb)%omx
                ap(i)%omy   = anb(k,nb)%omy
                ap(i)%omz   = anb(k,nb)%omz
             endif
          endif
          !    enddo
          nb=2
          !    do k=1,pmax_nb(nb)
          !      if ( (mslv_nb(k,nb).eq.idp) ) then
          ! neighbor nb is master of particle mslv_nb(p,nb)
          call binsearch(idp,mslv_nb(1:pmax_nb(nb),nb),pmax_nb(nb),found_mstr,k)
          if(found_mstr) then
             dist = sqrt( (rightbound-anb(k,nb)%x)**2. + (frontbound-anb(k,nb)%y)**2. )
             if ( dist .lt. (radius+offset) ) then
                if(count_slve_loc.eq. 0 ) i = i+1
                ap(i)%mslv  = -mslv_nb(k,nb) ! myid is slave of particle mslv_nb(k,nb)
                ap(i)%nb(nb) = 1             ! neighbor nb of myid is particle's master
                count_slve_loc = count_slve_loc + 1
                isperiodx = 0.
                isperiody = 0.
                if (anb(k,nb)%x.lt.0.) isperiodx =  1.
                if (anb(k,nb)%x.gt.lx) isperiodx = -1.
                if (anb(k,nb)%y.lt.0.) isperiody =  1.
                if (anb(k,nb)%y.gt.ly) isperiody = -1.
                ap(i)%x = anb(k,nb)%x+isperiodx*lx
                ap(i)%y = anb(k,nb)%y+isperiody*ly
                ap(i)%z = anb(k,nb)%z
                ap(i)%theta = anb(k,nb)%theta
                ap(i)%phi   = anb(k,nb)%phi
                ap(i)%u     = anb(k,nb)%u
                ap(i)%v     = anb(k,nb)%v
                ap(i)%w     = anb(k,nb)%w
                ap(i)%omx   = anb(k,nb)%omx
                ap(i)%omy   = anb(k,nb)%omy
                ap(i)%omz   = anb(k,nb)%omz
             endif
          endif
          !    enddo
          nb=3
          !    do k=1,pmax_nb(nb)
          !      if ( (mslv_nb(k,nb).eq.idp) ) then
          ! neighbor nb is master of particle mslv_nb(p,nb)
          call binsearch(idp,mslv_nb(1:pmax_nb(nb),nb),pmax_nb(nb),found_mstr,k)
          if(found_mstr) then
             dist = sqrt( (frontbound-anb(k,nb)%y)**2. )
             if ( dist .lt. (radius+offset) ) then
                if(count_slve_loc.eq. 0 ) i = i+1
                ap(i)%mslv  = -mslv_nb(k,nb) ! myid is slave of particle mslv_nb(k,nb)
                ap(i)%nb(nb) = 1             ! neighbor nb of myid is particle's master
                count_slve_loc = count_slve_loc + 1
                isperiodx = 0.
                isperiody = 0.
                if (anb(k,nb)%x.lt.0.) isperiodx =  1.
                if (anb(k,nb)%x.gt.lx) isperiodx = -1.
                if (anb(k,nb)%y.lt.0.) isperiody =  1.
                if (anb(k,nb)%y.gt.ly) isperiody = -1.
                ap(i)%x = anb(k,nb)%x+isperiodx*lx
                ap(i)%y = anb(k,nb)%y+isperiody*ly
                ap(i)%z = anb(k,nb)%z
                ap(i)%theta = anb(k,nb)%theta
                ap(i)%phi   = anb(k,nb)%phi
                ap(i)%u     = anb(k,nb)%u
                ap(i)%v     = anb(k,nb)%v
                ap(i)%w     = anb(k,nb)%w
                ap(i)%omx   = anb(k,nb)%omx
                ap(i)%omy   = anb(k,nb)%omy
                ap(i)%omz   = anb(k,nb)%omz
             endif
          endif
          !    enddo
          nb=4
          !    do k=1,pmax_nb(nb)
          !      if ( (mslv_nb(k,nb).eq.idp) ) then
          ! neighbor nb is master of particle mslv_nb(p,nb)
          call binsearch(idp,mslv_nb(1:pmax_nb(nb),nb),pmax_nb(nb),found_mstr,k)
          if(found_mstr) then
             dist = sqrt( (leftbound-anb(k,nb)%x)**2. + (frontbound-anb(k,nb)%y)**2. )
             if ( dist .lt. (radius+offset) ) then
                if(count_slve_loc.eq. 0 ) i = i+1
                ap(i)%mslv  = -mslv_nb(k,nb) ! myid is slave of particle mslv_nb(k,nb)
                ap(i)%nb(nb) = 1             ! neighbor nb of myid is particle's master
                count_slve_loc = count_slve_loc + 1
                isperiodx = 0.
                isperiody = 0.
                if (anb(k,nb)%x.lt.0.) isperiodx =  1.
                if (anb(k,nb)%x.gt.lx) isperiodx = -1.
                if (anb(k,nb)%y.lt.0.) isperiody =  1.
                if (anb(k,nb)%y.gt.ly) isperiody = -1.
                ap(i)%x = anb(k,nb)%x+isperiodx*lx
                ap(i)%y = anb(k,nb)%y+isperiody*ly
                ap(i)%z = anb(k,nb)%z
                ap(i)%theta = anb(k,nb)%theta
                ap(i)%phi   = anb(k,nb)%phi
                ap(i)%u     = anb(k,nb)%u
                ap(i)%v     = anb(k,nb)%v
                ap(i)%w     = anb(k,nb)%w
                ap(i)%omx   = anb(k,nb)%omx
                ap(i)%omy   = anb(k,nb)%omy
                ap(i)%omz   = anb(k,nb)%omz
             endif
          endif
          !    enddo
          if(count_slve_loc.ne.0) count_slve = count_slve + 1
       endif
    enddo
    !
    ! the new value of pmax is:
    !
    pmax = count_mstr + count_slve
    npmstr = count_mstr
    !
    !check if number of masters yield np
    !
    call MPI_ALLREDUCE(count_mstr,count_mstr_all,1,MPI_INTEGER,MPI_SUM,comm_cart,error)
    if(myid.eq.0) write(6,*) 'num of particles = num of masters?',np,' = ', count_mstr_all
    if (count_mstr_all.ne.np) then
       write(6,*) 'Fatal error in initialisation of particle positions!'
       write(6,*) 'Program aborted...'
       call mpi_abort(comm_cart,error,error)
       stop
    elseif(pmax.gt.npmax) then
       write(6,*) 'Size of local particle array of process ', myid, ' is too small!'
       write(6,*) 'Program aborted... (later I will simply write a warning and allocate more memory)'
       call mpi_abort(comm_cart,error,error)
       stop
    endif
    !
    !write new master/slave configuration to a file for debugging purposes
    !
    !write(rankpr,'(i3.3)') myid
    !open(25,file=datadir//'mslv'//rankpr//'.txt')
    !do p=1,pmax
    !  idp = abs(ap(p)%mslv)
    !  if (ap(p)%mslv .gt. 0) then
    !    counter = 0
    !    do i=1,8
    !      if (ap(p)%nb(i) .eq. 1) then
    !        counter = counter+1
    !        write(25,'(I4,A1,I5,A1,I5,A1,I2,A1,I2,2E16.8)') &
    !              myid,' ',idp,' ',ap(p)%mslv,' ',i,' ',neighbor(i),ap(p)%x,ap(p)%y
    !!        write(6,'(A29,I4,A1,I5,A1,I5,A1,I2,A1,I2,2E16.8)') 'rank,p,pms,nbr,ranknbr,x,y = ', &
    !!              myid,' ',idp,' ',ap(p)%mslv,' ',i,' ',neighbor(i),ap(p)%x,ap(p)%y
    !      endif
    !    enddo
    !    !in case of no overlap with any neighbor
    !    if (counter .eq. 0) then
    !      write(25,'(I4,A1,I5,A1,I5,A1,I2,A1,I2,2E16.8)') &
    !            myid,' ',idp,' ',ap(p)%mslv,' ',99,' ',99,ap(p)%x,ap(p)%y
    !!      write(6,'(A29,I4,A1,I5,A1,I5,A1,I2,A1,I2,2E16.8)') 'rank,p,pms,nbr,ranknbr,x,y = ', &
    !!              myid,' ',idp,' ',ap(p)%mslv,' ',99,' ',99,ap(p)%x,ap(p)%y
    !    endif
    !  endif
    !  if (ap(p)%mslv .lt. 0) then
    !    do i=1,8
    !      if (ap(p)%nb(i) .eq. 1) then
    !        write(25,'(I4,A1,I5,A1,I5,A1,I2,A1,I2,2E16.8)') &
    !              myid,' ',idp,' ',ap(p)%mslv,' ',i,' ',neighbor(i),ap(p)%x,ap(p)%y
    !!        write(6,'(A29,I4,A1,I5,A1,I5,A1,I2,A1,I2,2E16.8)') 'rank,p,pms,nbr,ranknbr,x,y = ', &
    !!              myid,' ',idp,' ',ap(p)%mslv,' ',i,' ',neighbor(i),ap(p)%x,ap(p)%y
    !      endif
    !    enddo
    !  endif
    !enddo
    !close(25)
    !
    ! masters: new positions and velocities of Lagrangian forcing points
    !
    !$omp workshare
    nla(:) = 0 ! set to zero from 1 to npmax
    !$omp end workshare
!
    !$omp parallel default(none)               &
    !$omp&shared(ap,nla,pmax)                  &
    !$omp&shared(radfp,phirc,thetarc)          &
    !$omp&shared(boundleftmyid,boundfrontmyid) &
    !$omp&private(p,l,ll,coorxfp,cooryfp,xfploc,yfploc,zfploc)    &
    !$omp&private(isperiodx,isperiody,isout)
    !$omp do 
    do p=1,pmax
       if (ap(p)%mslv .ne. 0) then
          ! myid is master of particle ap(p)%mslv
          ll = 0
          do l=1,NL
             xfploc = ap(p)%x + radfp*sin(ap(p)%theta*0.+thetarc(l))*cos(ap(p)%phi*0.+phirc(l))
             yfploc = ap(p)%y + radfp*sin(ap(p)%theta*0.+thetarc(l))*sin(ap(p)%phi*0.+phirc(l))
             zfploc = ap(p)%z + radfp*cos(ap(p)%theta*0.+thetarc(l))
             isperiodx = 0.
             isperiody = 0.
             if (xfploc.lt.0.+0.5*dx) isperiodx =  1.
             if (xfploc.ge.lx+0.5*dx) isperiodx = -1.
             if (yfploc.lt.0.+0.5*dy) isperiody =  1.
             if (yfploc.ge.ly+0.5*dy) isperiody = -1.
             isout = .false.
             coorxfp = (xfploc+isperiodx*lx-boundleftmyid )*dxi
             if( nint(coorxfp).lt.1 .or. nint(coorxfp) .gt.imax ) isout = .true.
             cooryfp = (yfploc+isperiody*ly-boundfrontmyid)*dyi
             if( nint(cooryfp).lt.1 .or. nint(cooryfp) .gt.jmax ) isout = .true.
!
             if(.not.isout) then
               ll = ll + 1
               ap(p)%xfp(ll) = xfploc
               ap(p)%yfp(ll) = yfploc
               ap(p)%zfp(ll) = zfploc
               ap(p)%ul(ll) = ap(p)%u + ap(p)%omy*(ap(p)%zfp(ll)-ap(p)%z) &
                                      - ap(p)%omz*(ap(p)%yfp(ll)-ap(p)%y)
               ap(p)%vl(ll) = ap(p)%v + ap(p)%omz*(ap(p)%xfp(ll)-ap(p)%x) &
                                      - ap(p)%omx*(ap(p)%zfp(ll)-ap(p)%z)
               ap(p)%wl(ll) = ap(p)%w + ap(p)%omx*(ap(p)%yfp(ll)-ap(p)%y) &
                                      - ap(p)%omy*(ap(p)%xfp(ll)-ap(p)%x)
             endif
          enddo
          nla(p) = ll
       endif
    enddo
    !$omp end parallel
!
!    do p=1,pmax
!       !  if (ap(p)%mslv .gt. 0) then
!       if (ap(p)%mslv .ne. 0) then
!          do l=1,NL
!             ap(p)%xfp(l) = ap(p)%x + radfp*sin(ap(p)%theta*0.+thetarc(l))*cos(ap(p)%phi*0.+phirc(l))
!             ap(p)%yfp(l) = ap(p)%y + radfp*sin(ap(p)%theta*0.+thetarc(l))*sin(ap(p)%phi*0.+phirc(l))
!             ap(p)%zfp(l) = ap(p)%z + radfp*cos(ap(p)%theta*0.+thetarc(l))
!             ap(p)%ul(l) = ap(p)%u + ap(p)%omy*(ap(p)%zfp(l)-ap(p)%z) &
!                  - ap(p)%omz*(ap(p)%yfp(l)-ap(p)%y)
!             ap(p)%vl(l) = ap(p)%v + ap(p)%omz*(ap(p)%xfp(l)-ap(p)%x) &
!                  - ap(p)%omx*(ap(p)%zfp(l)-ap(p)%z)
!             ap(p)%wl(l) = ap(p)%w + ap(p)%omx*(ap(p)%yfp(l)-ap(p)%y) &
!                  - ap(p)%omy*(ap(p)%xfp(l)-ap(p)%x)
!          enddo
!          !  else
!          !    do l=1,NL
!          !      ap(p)%xfp(l) = 0.
!          !      ap(p)%yfp(l) = 0.
!          !      ap(p)%zfp(l) = 0.
!          !      ap(p)%ul(l)  = 0.
!          !      ap(p)%vl(l)  = 0.
!          !      ap(p)%wl(l)  = 0.
!          !    enddo
!       endif
!    enddo
    !
    return
  end subroutine intgr_nwtn_eulr
  !

!  subroutine binsearch(ival,array,idim,found,index)
!    integer, intent(in), dimension(1:) :: array
!    integer, intent(in) :: ival,idim
!    logical, intent(out) :: found
!    integer, intent(out) :: index
!    integer :: start,finish,range,mid

!    start = 1
!    finish = idim
!    range = finish-start
!    mid = (start+finish)/2
!    do while( abs(array(mid)) .ne. ival .and. range .gt.  0)
!       if (ival .gt. abs(array(mid))) then
!          start = mid + 1
!       else
!          finish = mid - 1
!       endif
!       range = finish - start
!       mid = (start + finish)/2
!    enddo
!    if(array(mid).ne.ival) then
!       found = .false.
!       index = 0
!    else
!       found = .true.
!       index = mid
!    endif
!    return
!  end subroutine binsearch


  subroutine binsearch(key,array,idim,found,index) 
    !! this is not a real binary search !!
    ! it searches key with absolute value, then checks its positivity
    integer, intent(in), dimension(1:) :: array
    integer, intent(in) :: key,idim
    logical, intent(out) :: found
    integer, intent(out) :: index
    integer imin, imax, imid

    found = .false.
    index = 0
    imin = 1
    imax = idim    

    ! continue searching while [imin,imax] is not empty
    do while (imin .le. imax)
       !calculate the midpoint for roughly equal partition
       imid = (imin + imax)/2;
       if(abs(array(imid)) .eq. key) then
          !key found at index imid 
          if (array(imid) .eq. key) then
             !! this is what makes it not a real binary search !!
             index = imid
             found = .true.
          endif
          return
          ! determine which subarray to search
       else if (abs(array(imid)) .lt. key) then
          !change min index to search upper subarray
          imin = imid + 1;
       else          
          !change max index to search lower subarray
          imax = imid - 1;
       end if
    end do
    return
  end subroutine binsearch

  subroutine sumrk3
    implicit none
    integer :: p
    !!
    !! accumulate values over RK3 substeps
    !!
!$omp parallel default(shared) &
!$omp&private(p)
!$omp do schedule(dynamic)
    do p=1,pmax
      !
      !initialisation at first RK sub-step
      !
      if (rkiter .eq. 1) then
        if (ap(p)%mslv .gt. 0) then
          rkp(p)%fxltot    = 0. 
          rkp(p)%fyltot    = 0. 
          rkp(p)%fzltot    = 0. 
          rkp(p)%torqxltot = 0. 
          rkp(p)%torqyltot = 0. 
          rkp(p)%torqzltot = 0. 
          rkp(p)%torqtheta = 0. 
          rkp(p)%colfx     = 0. 
          rkp(p)%colfy     = 0. 
          rkp(p)%colfz     = 0. 
          rkp(p)%coltx     = 0. 
          rkp(p)%colty     = 0. 
          rkp(p)%coltz     = 0. 
          rkp(p)%dudt      = 0. 
          rkp(p)%dvdt      = 0. 
          rkp(p)%dwdt      = 0. 
          rkp(p)%domxdt   = 0.
          rkp(p)%domydt   = 0.
          rkp(p)%domzdt   = 0.
        endif
      endif
      !
      ! accumulation
      !
      if (ap(p)%mslv .gt. 0) then
        rkp(p)%fxltot = rkp(p)%fxltot + &
            (-1.)*ap(p)%fxltot + (ap(p)%intu-op(p)%intu)/dt  ! multiplied with dVlagr!
        rkp(p)%fyltot = rkp(p)%fyltot + &
            (-1.)*ap(p)%fyltot + (ap(p)%intv-op(p)%intv)/dt  ! multiplied with dVlagr!
        rkp(p)%fzltot = rkp(p)%fzltot + &
            (-1.)*ap(p)%fzltot + (ap(p)%intw-op(p)%intw)/dt  ! multiplied with dVlagr!
        rkp(p)%torqxltot         = rkp(p)%torqxltot     + &
            (-1.)*ap(p)%torqxltot + (ap(p)%intomx-op(p)%intomx)/dt
        rkp(p)%torqyltot         = rkp(p)%torqyltot     + &
            (-1.)*ap(p)%torqyltot + (ap(p)%intomy-op(p)%intomy)/dt
        rkp(p)%torqzltot         = rkp(p)%torqzltot     + &
            (-1.)*ap(p)%torqzltot + (ap(p)%intomz-op(p)%intomz)/dt
        rkp(p)%torqtheta     = rkp(p)%torqtheta + & 
            (-1.)*(-dpdy_new*dVlagr)!(-1.)*ap(p)%torqtheta ! note: not corrected for inertia of fluid inside sphere
        rkp(p)%colfx = rkp(p)%colfx   + &
                             0.5*rkcoeffab(rkiter)*(ap(p)%colfx+op(p)%colfx)
        rkp(p)%colfy = rkp(p)%colfy   + &
                             0.5*rkcoeffab(rkiter)*(ap(p)%colfy+op(p)%colfy)
        rkp(p)%colfz = rkp(p)%colfz   + &
                             0.5*rkcoeffab(rkiter)*(ap(p)%colfz+op(p)%colfz)
        rkp(p)%coltx = rkp(p)%coltx   + &
                             0.5*rkcoeffab(rkiter)*(ap(p)%coltx+op(p)%coltx)
        rkp(p)%colty = rkp(p)%colty   + &
                             0.5*rkcoeffab(rkiter)*(ap(p)%colty+op(p)%colty)
        rkp(p)%coltz = rkp(p)%coltz   + &
                             0.5*rkcoeffab(rkiter)*(ap(p)%coltz+op(p)%coltz)
        rkp(p)%dudt  = rkp(p)%dudt + (ap(p)%u-op(p)%u)/dt
        rkp(p)%dvdt  = rkp(p)%dvdt + (ap(p)%v-op(p)%v)/dt
        rkp(p)%dwdt  = rkp(p)%dwdt + (ap(p)%w-op(p)%w)/dt
        ! check acceleration calculation ! (SS)
        !   ! Only compute acceleration at the last RK substep (rkiter=3)
      !   ! to get meaningful physical acceleration based on full timestep
      !   if (rkiter .eq. 3) then
      !     rkp(p)%dudt  = (ap(p)%u-op(p)%u)/dt  ! total velocity change over full timestep
      !     rkp(p)%dvdt  = (ap(p)%v-op(p)%v)/dt
      !     rkp(p)%dwdt  = (ap(p)%w-op(p)%w)/dt
      !   endif

      endif
    enddo
    !$omp end parallel
    !
    return
  end subroutine sumrk3
  subroutine transfer_sumrk3
      implicit none
      integer :: p
!$omp parallel default(shared) &
!$omp&private(p)
!$omp do schedule(dynamic)
    do p=1,pmax
      srkp(p)%fxltot    = rkp(p)%fxltot
      srkp(p)%fyltot    = rkp(p)%fyltot
      srkp(p)%fzltot    = rkp(p)%fzltot
      srkp(p)%torqxltot = rkp(p)%torqxltot
      srkp(p)%torqyltot = rkp(p)%torqyltot
      srkp(p)%torqzltot = rkp(p)%torqzltot
      srkp(p)%torqtheta = rkp(p)%torqtheta
      srkp(p)%colfx     = rkp(p)%colfx
      srkp(p)%colfy     = rkp(p)%colfy
      srkp(p)%colfz     = rkp(p)%colfz
      srkp(p)%coltx     = rkp(p)%coltx
      srkp(p)%colty     = rkp(p)%colty
      srkp(p)%coltz     = rkp(p)%coltz
      srkp(p)%dudt      = rkp(p)%dudt
      srkp(p)%dvdt      = rkp(p)%dvdt
      srkp(p)%dwdt      = rkp(p)%dwdt
      srkp(p)%domxdt    = rkp(p)%domxdt
      srkp(p)%domydt    = rkp(p)%domydt
      srkp(p)%domzdt    = rkp(p)%domzdt
    enddo
    do p=1,npmax
      rkp(p)%fxltot    = 0. 
      rkp(p)%fyltot    = 0. 
      rkp(p)%fzltot    = 0. 
      rkp(p)%torqxltot = 0. 
      rkp(p)%torqyltot = 0. 
      rkp(p)%torqzltot = 0. 
      rkp(p)%torqtheta = 0. 
      rkp(p)%colfx     = 0. 
      rkp(p)%colfy     = 0. 
      rkp(p)%colfz     = 0. 
      rkp(p)%coltx     = 0. 
      rkp(p)%colty     = 0. 
      rkp(p)%coltz     = 0. 
      rkp(p)%dudt      = 0. 
      rkp(p)%dvdt      = 0. 
      rkp(p)%dwdt      = 0. 
      rkp(p)%domxdt   = 0.
      rkp(p)%domydt   = 0.
      rkp(p)%domzdt   = 0.
    enddo
    !$omp end parallel
  end subroutine transfer_sumrk3
  subroutine sumrk3_newmaster(q,p)
      implicit none
    integer :: p,q
    rkp(q)%fxltot    = srkp(p)%fxltot
    rkp(q)%fyltot    = srkp(p)%fyltot
    rkp(q)%fzltot    = srkp(p)%fzltot
    rkp(q)%torqxltot = srkp(p)%torqxltot
    rkp(q)%torqyltot = srkp(p)%torqyltot
    rkp(q)%torqzltot = srkp(p)%torqzltot
    rkp(q)%torqtheta = srkp(p)%torqtheta
    rkp(q)%colfx     = srkp(p)%colfx
    rkp(q)%colfy     = srkp(p)%colfy
    rkp(q)%colfz     = srkp(p)%colfz
    rkp(q)%coltx     = srkp(p)%coltx
    rkp(q)%colty     = srkp(p)%colty
    rkp(q)%coltz     = srkp(p)%coltz
    rkp(q)%dudt      = srkp(p)%dudt
    rkp(q)%dvdt      = srkp(p)%dvdt
    rkp(q)%dwdt      = srkp(p)%dwdt
  end subroutine sumrk3_newmaster
  !
end module mod_intgr_nwtn_eulr
