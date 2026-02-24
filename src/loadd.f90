module mod_loadflds
  use decomp_2d
  use decomp_2d_io
  use mod_param 
  use mod_common
  use mod_common_mpi
  implicit none
  private
  public loadflds,amplify_vel_fluctuations
contains
  !
  subroutine loadflds(in,nr)
    implicit none
    integer :: in,nr
    integer :: fh
    integer(kind=MPI_OFFSET_KIND) :: filesize,disp
    character(len=7) :: istepchar
    real, dimension(3) :: fldinfo
    real, dimension(imax,jmax,kmax) :: temp ! should be changed to a pointer
    !
    if (in.eq.0) then
       write(istepchar,'(i7.7)') nr!0
       call MPI_FILE_OPEN(MPI_COMM_WORLD, 'fld'//istepchar, &
            MPI_MODE_RDONLY, MPI_INFO_NULL,fh, error)
       disp = 0_MPI_OFFSET_KIND
       call decomp_2d_read_var(fh,disp,3,temp) ! read fields
       unew(1:imax,1:jmax,1:kmax) = temp
       call decomp_2d_read_var(fh,disp,3,temp)
       vnew(1:imax,1:jmax,1:kmax) = temp
       call decomp_2d_read_var(fh,disp,3,temp)
       wnew(1:imax,1:jmax,1:kmax) = temp
       call decomp_2d_read_var(fh,disp,3,temp)
       pnew(1:imax,1:jmax,1:kmax) = temp
       call decomp_2d_read_var(fh,disp,3,temp)
       cnew(1:imax,1:jmax,1:kmax) = temp
       
       call decomp_2d_read_scalar(fh,disp,3,fldinfo) ! read scalar
       time = fldinfo(1)
       nr = int(fldinfo(2))
       dt = fldinfo(3)
       call MPI_FILE_CLOSE(fh,error)
    endif
    !
    if (in.eq.1) then
       write(istepchar,'(i7.7)') nr !0
       fldinfo = (/time,1.*nr,dt/)
       call MPI_FILE_OPEN(MPI_COMM_WORLD, 'fld'//istepchar, &
            MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, error)
       filesize = 0_MPI_OFFSET_KIND
       call MPI_FILE_SET_SIZE(fh,filesize,error)  ! guarantee overwriting
       disp = 0_MPI_OFFSET_KIND
       temp = unew(1:imax,1:jmax,1:kmax)
       call decomp_2d_write_var(fh,disp,3,temp)
       temp = vnew(1:imax,1:jmax,1:kmax)
       call decomp_2d_write_var(fh,disp,3,temp)
       temp = wnew(1:imax,1:jmax,1:kmax)
       call decomp_2d_write_var(fh,disp,3,temp)
       temp = pnew(1:imax,1:jmax,1:kmax)
       call decomp_2d_write_var(fh,disp,3,temp)
       temp = cnew(1:imax,1:jmax,1:kmax)

       call decomp_2d_write_var(fh,disp,3,temp)
       call decomp_2d_write_scalar(fh,disp,3,fldinfo)
       call MPI_FILE_CLOSE(fh,error)
    endif
    !
    return
  end subroutine loadflds
  !
  subroutine amplify_vel_fluctuations(factor)
    implicit none
    real, dimension(kmax) :: um,vm,wm,sum_all
    real, intent(in) :: factor
    integer :: i,j,k
    !
    do k=1,kmax
       um(k) = 0.
       vm(k) = 0.
       wm(k) = 0.
       do j=1,jmax
          do i=1,imax
             um(k)  = um(k)  + unew(i,j,k)
             vm(k)  = vm(k)  + vnew(i,j,k)
             wm(k)  = wm(k)  + wnew(i,j,k)
          enddo
       enddo
    enddo
    !
    call mpi_allreduce(um(1),sum_all(1),kmax,mpi_real8,mpi_sum,comm_cart,error)
    um(:) = sum_all(:)/(1.*itot*jtot)
    call mpi_allreduce(vm(1),sum_all(1),kmax,mpi_real8,mpi_sum,comm_cart,error)
    vm(:) = sum_all(:)/(1.*itot*jtot)
    call mpi_allreduce(wm(1),sum_all(1),kmax,mpi_real8,mpi_sum,comm_cart,error)
    wm(:) = sum_all(:)/(1.*itot*jtot)
    !
    do k=1,kmax
       do j=1,jmax
          do i=1,imax
             unew(i,j,k) = (unew(i,j,k) - um(k))*factor
             vnew(i,j,k) = (vnew(i,j,k) - vm(k))*factor + vm(k)
             wnew(i,j,k) = (wnew(i,j,k) - wm(k))*factor
          enddo
       enddo
    enddo
    !
    return
  end subroutine amplify_vel_fluctuations
  !
end module mod_loadflds
!
module mod_loadpart
  use mpi
  use mod_param
  use mod_common
  use mod_common_mpi
  implicit none
  private
  public loadpart
contains
  !
  subroutine loadpart(in,nr)
    implicit none
    type particle_restart
       real :: x,y,z,theta,phi, &
            u,v,w, &
            omx,omy,omz,omtheta, &
            intu,intv,intw, &
            intomx,intomy,intomz, &
            flag1,flag2, &
            colfx,colfy,colfz, &
            coltx,colty,coltz, &
            fxltot, fyltot, fzltot !
       real, dimension(nqmax) :: dx,dy,dz, &
            dxt,dyt,dzt, &
            dut,dvt,dwt, &
            psi ! 10*nqmax
       ! total ammount of reals to be communicated: 24+10*nqmax
       real :: qmax,idp ! 2
       real, dimension(nqmax) :: firstc ! 1*nqmax
    end type particle_restart
    type(particle_restart), allocatable, dimension(:) :: glob
    integer, parameter :: skipr = 27+10*nqmax, skipi = 2+nqmax, & 
         skip = skipr+skipi
    integer,dimension(0:dims(1)*dims(2)-1) :: npmstr_glob,npmstr_glob_all
    ! skipi is is the same as in the common flie + the global id of the particle (COBP)
    integer i,j,p,idp
    integer :: fh
    integer     in,nr
    integer :: proccoords(1:ndims),procrank
    real :: leftbound,rightbound,frontbound,backbound
    real :: xdum
    integer :: lenr
    real :: dist,angle
    real :: xp,yp
    real :: ax
    real :: ay
    integer :: count_mstr_all
    character(len=7) :: istepchar
    character(len=4) rankpr
    character(len=11) :: tempstr2
    integer :: counter,count_slve_loc
    integer(kind=MPI_OFFSET_KIND) :: filesize,disp
    logical :: found_mstr
    integer, dimension(send_int) :: itemp
    integer :: mydisp
    integer, dimension(np) :: rkmstr_glob,rkmstr_fake,locid_fake,locid_true
    integer, dimension(np) :: rkmstr_glob_all,rkmstr_fake_all,locid_fake_all,locid_true_all
    real, dimension(np) :: xcglob,ycglob,zcglob
    real, dimension(np) :: xcglob_all,ycglob_all,zcglob_all
    integer :: nrrequests
    integer :: arrayrequests(1:3)
    integer :: arraystatuses(MPI_STATUS_SIZE,1:3)
    integer :: tag,nb,nbsend,nbrecv
    integer :: premain
    !
    !inquire (iolength=lenr) xdum
    lenr = sizeof(xdum)
    write(istepchar,'(i7.7)') nr !0
    !
    if (in.eq.0) then
       pmax    = np/product(dims)
       premain = np - pmax*product(dims)
       if( myid.lt.premain ) pmax = pmax + 1
       allocate(glob(pmax))
       call MPI_FILE_OPEN(MPI_COMM_WORLD, 'allpartdata'//istepchar, &
            MPI_MODE_RDONLY, MPI_INFO_NULL,fh, error)
       npmstr_glob(:) = 0
       npmstr_glob(myid) = pmax
       call MPI_ALLREDUCE(npmstr_glob(0),npmstr_glob_all(0),product(dims),MPI_INTEGER,MPI_SUM,comm_cart,error)
       mydisp = 0
       if(myid.ne.0) mydisp = sum(npmstr_glob_all(0:myid-1))
       disp = mydisp*skip*lenr
       call MPI_FILE_SET_VIEW(fh, disp, MPI_REAL8,MPI_REAL8, 'native', & 
            MPI_INFO_NULL, error)
       call MPI_FILE_READ(fh,glob(1)%x,skip*pmax,MPI_REAL8,MPI_STATUS_IGNORE,error)
       call MPI_FILE_CLOSE(fh,error)
       rkmstr_fake(:) = 0
       locid_fake(:) = 0
       locid_true(:) = 0
       rkmstr_glob(:) = 0
       xcglob(:) = 0.
       ycglob(:) = 0.
       zcglob(:) = 0.
       do p=1,pmax
          idp = nint(glob(p)%idp)
          rkmstr_fake(idp) = myid
          locid_fake(idp) = p
          if (glob(p)%x.lt.0..or.glob(p)%x.gt.lx .or. &
               glob(p)%y.lt.0..or.glob(p)%y.gt.ly .or. &
               glob(p)%z.lt.0..or.glob(p)%z.gt.lz) then
             if (myid.eq.0) then
                write(6,*) 'Fatal error in initialisation of particle positions - '
                write(6,*) 'particle outside the domain!'
                write(6,*) 'Program aborted...'
             endif
             call mpi_finalize(error)
             stop
          endif
          ax = 0.5
          ay = 0.5
          if (glob(p)%x.eq.lx) ax = 0.51
          if (glob(p)%x.eq.0.) ax = 0.49
          if (glob(p)%y.eq.ly) ay = 0.51
          if (glob(p)%y.eq.0.) ay = 0.49
          proccoords(1) = nint(glob(p)%x*dims(1)/lx - ax)
          proccoords(2) = nint(glob(p)%y*dims(2)/ly - ay)
          call MPI_CART_RANK(comm_cart,proccoords,procrank,error)
          rkmstr_glob(idp) = procrank
          xcglob(idp) = glob(p)%x
          ycglob(idp) = glob(p)%y
          zcglob(idp) = glob(p)%z
       enddo
       call MPI_ALLREDUCE(rkmstr_glob(1),rkmstr_glob_all(1),np,MPI_INTEGER,MPI_SUM,comm_cart,error)
       call MPI_ALLREDUCE(rkmstr_fake(1),rkmstr_fake_all(1),np,MPI_INTEGER,MPI_SUM,comm_cart,error)
       call MPI_ALLREDUCE(locid_fake(1) ,locid_fake_all(1) ,np,MPI_INTEGER,MPI_SUM,comm_cart,error)
       call MPI_ALLREDUCE(xcglob(1) ,xcglob_all(1) ,np,MPI_REAL8,MPI_SUM,comm_cart,error)
       call MPI_ALLREDUCE(ycglob(1) ,ycglob_all(1) ,np,MPI_REAL8,MPI_SUM,comm_cart,error)
       call MPI_ALLREDUCE(zcglob(1) ,zcglob_all(1) ,np,MPI_REAL8,MPI_SUM,comm_cart,error)
       rkmstr_glob(:) = rkmstr_glob_all(:)
       rkmstr_fake(:) = rkmstr_fake_all(:)
       locid_fake(:)  = locid_fake_all(:)
       xcglob(:) = xcglob_all(:)
       ycglob(:) = ycglob_all(:)
       zcglob(:) = zcglob_all(:)
       i = 0
       do idp = 1,np
          nrrequests = 0
          if(rkmstr_glob(idp).eq.myid) then
             i = i + 1
             nrrequests = nrrequests + 1
             CALL MPI_IRECV(sp(i)%x,send_real,MPI_REAL8,rkmstr_fake(idp),idp,comm_cart,arrayrequests((nrrequests-1)*2+1),error)
             CALL MPI_IRECV(sp(i)%qmax,send_int,MPI_INTEGER,rkmstr_fake(idp),idp+np,comm_cart,arrayrequests((nrrequests-1)*2+2),error)
             sp(i)%mslv = idp
             locid_true(idp) = i
          endif
          if(rkmstr_fake(idp).eq.myid) then
             p = locid_fake(idp)
             nrrequests = nrrequests + 1
             CALL MPI_ISEND(glob(p)%x,send_real,MPI_REAL8,rkmstr_glob(idp),idp,comm_cart,arrayrequests((nrrequests-1)*2+1),error)
             itemp(1) = nint(glob(p)%qmax)
             itemp(2:nqmax+1) =  nint(glob(p)%firstc(:))
             CALL MPI_ISEND(itemp(1),send_int,MPI_INTEGER,rkmstr_glob(idp),idp+np,comm_cart,arrayrequests((nrrequests-1)*2+2),error)
          endif
          nrrequests = nrrequests*2
          call MPI_WAITALL(nrrequests,arrayrequests,arraystatuses,error)
       enddo
       call MPI_ALLREDUCE(locid_true(1) ,locid_true_all(1) ,np,MPI_INTEGER,MPI_SUM,comm_cart,error)
       locid_true(:)  = locid_true_all(:)
       pmax = i
       ! 
       ! Determine master and slave processes for each particle.
       !
       ! initialisation
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
       ap(1:npmax)%qmax = 0
       forall (p=1:npmax)
          ap(p)%dx(1:nqmax) = 0. 
          ap(p)%dy(1:nqmax) = 0.
          ap(p)%dz(1:nqmax) = 0.
          ap(p)%dxt(1:nqmax) = 0.
          ap(p)%dyt(1:nqmax) = 0.
          ap(p)%dzt(1:nqmax) = 0.
          ap(p)%firstc(1:nqmax) = 0
          ap(p)%nb(1:8) = 0
       end forall
       i = 0
       npmstr = 0
       do idp=1,np
          ax = 0.5
          ay = 0.5
          if (xcglob(idp).eq.lx) ax = 0.51
          if (xcglob(idp).eq.0.) ax = 0.49
          if (ycglob(idp).eq.ly) ay = 0.51
          if (ycglob(idp).eq.0.) ay = 0.49
          proccoords(1) = nint(xcglob(idp)*dims(1)/lx - ax)
          proccoords(2) = nint(ycglob(idp)*dims(2)/ly - ay)
          leftbound     = (proccoords(1)  )*lx/(1.*dims(1)) ! left  boundary of particle's master
          rightbound    = (proccoords(1)+1)*lx/(1.*dims(1)) ! right boundary of particle's master
          frontbound    = (proccoords(2)  )*ly/(1.*dims(2)) ! front boundary of particle's master
          backbound     = (proccoords(2)+1)*ly/(1.*dims(2)) ! back  boundary of particle's master
          if(rkmstr_glob(idp).eq.myid) then
             npmstr = npmstr + 1
             i = i + 1
             p = locid_true(idp)
             ap(i)%mslv      = sp(p)%mslv
             ap(i)%x         = sp(p)%x
             ap(i)%y         = sp(p)%y
             ap(i)%z         = sp(p)%z
             ap(i)%theta     = sp(p)%theta
             ap(i)%phi       = sp(p)%phi
             ap(i)%u         = sp(p)%u
             ap(i)%v         = sp(p)%v
             ap(i)%w         = sp(p)%w
             ap(i)%omx       = sp(p)%omx
             ap(i)%omy       = sp(p)%omy
             ap(i)%omz       = sp(p)%omz
             ap(i)%omtheta   = sp(p)%omtheta
             ap(i)%intu      = sp(p)%intu
             ap(i)%intv      = sp(p)%intv
             ap(i)%intw      = sp(p)%intw
             ap(i)%intomx    = sp(p)%intomx
             ap(i)%intomy    = sp(p)%intomy
             ap(i)%intomz    = sp(p)%intomz
             ap(i)%colfx     = sp(p)%colfx
             ap(i)%colfy     = sp(p)%colfy
             ap(i)%colfz     = sp(p)%colfz
             ap(i)%coltx     = sp(p)%coltx
             ap(i)%colty     = sp(p)%colty
             ap(i)%coltz     = sp(p)%coltz
             ap(i)%dx(1:nqmax) = sp(p)%dx(1:nqmax)
             ap(i)%dy(1:nqmax) = sp(p)%dy(1:nqmax)
             ap(i)%dz(1:nqmax) = sp(p)%dz(1:nqmax)
             ap(i)%dxt(1:nqmax) = sp(p)%dxt(1:nqmax)
             ap(i)%dyt(1:nqmax) = sp(p)%dyt(1:nqmax) 
             ap(i)%dzt(1:nqmax) = sp(p)%dzt(1:nqmax)
             ap(i)%qmax = sp(p)%qmax
             ap(i)%firstc(1:nqmax) = sp(p)%firstc(1:nqmax)
             ! neighbor 1
             if ( ap(i)%x .gt. (rightbound-(radius+offset)) ) then 
                ap(i)%nb(1) = 1 ! neighbor 1 is slave of particle ap(p)%mslv 
             endif
             ! neighbor 2
             dist = sqrt( (rightbound-ap(i)%x)**2. + (frontbound-ap(i)%y)**2. )
             if ( abs(dist) .lt. (radius+offset) ) then
                ap(i)%nb(2) = 1 ! neighbor 2 is slave of particle ap(p)%mslv
             endif
             ! neighbor 3
             if ( ap(i)%y .lt. (frontbound+(radius+offset)) ) then
                ap(i)%nb(3) = 1 ! neighbor 3 is slave of particle ap(p)%mslv
             endif
             ! neighbor 4
             dist = sqrt( (leftbound-ap(i)%x)**2. + (frontbound-ap(i)%y)**2. ) 
             if ( abs(dist) .lt. (radius+offset) ) then
                ap(i)%nb(4) = 1 ! neighbor 4 is slave of particle ap(p)%mslv
             endif
             ! neighbor 5
             if ( ap(i)%x .lt. (leftbound+(radius+offset)) ) then
                ap(i)%nb(5) = 1 ! neighbor 5 is slave of particle ap(p)%mslv
             endif
             ! neighbor 6
             dist = sqrt( (leftbound-ap(i)%x)**2. + (backbound-ap(i)%y)**2. )
             if ( abs(dist) .lt. (radius+offset) ) then
                ap(i)%nb(6) = 1 ! neighbor 6 is slave of particle ap(p)%mslv
             endif
             ! neighbor 7
             if ( ap(i)%y .gt. (backbound-(radius+offset)) ) then
                ap(i)%nb(7) = 1 ! neighbor 7 is slave of particle ap(p)%mslv
             endif
             ! neighbor 8
             dist = sqrt( (rightbound-ap(i)%x)**2. + (backbound-ap(i)%y)**2. )
             if ( abs(dist) .lt. (radius+offset) ) then
                ap(i)%nb(8) = 1 ! neighbor 8 is slave of particle ap(p)%mslv
             endif
          else
             count_slve_loc = 0
             ! neighbor 1 of particle's master
             proccoords(1) = nint( dims(1)*xcglob(idp)/lx - ax ) + 1
             proccoords(2) = nint( dims(2)*ycglob(idp)/ly - ay ) 
             call MPI_CART_RANK(comm_cart,proccoords,procrank,error)
             if (myid .eq. procrank) then
                if ( xcglob(idp) .gt. (rightbound-(radius+offset)) ) then
                   if(count_slve_loc.eq. 0 ) i = i+1
                   count_slve_loc = count_slve_loc + 1
                   ap(i)%mslv = -idp     ! myid is slave of particle abs(ap(p)%mslv)
                   ap(i)%nb(5) = 1      ! neighbor 5 of myid is particle's master
                endif
             endif
             ! neighbor 2 of particle's master
             proccoords(1) = nint( dims(1)*xcglob(idp)/lx - ax ) + 1
             proccoords(2) = nint( dims(2)*ycglob(idp)/ly - ay ) - 1
             call MPI_CART_RANK(comm_cart,proccoords,procrank,error)
             if (myid .eq. procrank) then
                dist = sqrt( (rightbound-xcglob(idp))**2. + (frontbound-ycglob(idp))**2. ) 
                if ( abs(dist) .lt. (radius+offset) ) then
                   if(count_slve_loc.eq. 0 ) i = i+1
                   count_slve_loc = count_slve_loc + 1
                   ap(i)%mslv = -idp     ! myid is slave of particle abs(ap(p)%mslv)
                   ap(i)%nb(6) = 1      ! neighbor 6 of myid is particle's master
                endif
             endif
             ! neighbor 3 of particle's master
             proccoords(1) = nint( dims(1)*xcglob(idp)/lx - ax )
             proccoords(2) = nint( dims(2)*ycglob(idp)/ly - ay ) - 1
             call MPI_CART_RANK(comm_cart,proccoords,procrank,error)
             if (myid .eq. procrank) then
                if ( ycglob(idp) .lt. (frontbound+(radius+offset)) ) then
                   if(count_slve_loc.eq. 0 ) i = i+1
                   count_slve_loc = count_slve_loc + 1
                   ap(i)%mslv = -idp     ! myid is slave of particle abs(ap(p)%mslv)
                   ap(i)%nb(7) = 1      ! neighbor 7 of myid is particle's master
                endif
             endif
             ! neighbor 4 of particle's master
             proccoords(1) = nint( dims(1)*xcglob(idp)/lx - ax ) - 1
             proccoords(2) = nint( dims(2)*ycglob(idp)/ly - ay ) - 1
             call MPI_CART_RANK(comm_cart,proccoords,procrank,error)
             if (myid .eq. procrank) then
                dist = sqrt( (leftbound-xcglob(idp))**2. + (frontbound-ycglob(idp))**2. )
                if ( abs(dist) .lt. (radius+offset) ) then
                   if(count_slve_loc.eq. 0 ) i = i+1
                   count_slve_loc = count_slve_loc + 1
                   ap(i)%mslv = -idp     ! myid is slave of particle abs(ap(p)%mslv)
                   ap(i)%nb(8) = 1      ! neighbor 8 of myid is particle's master
                endif
             endif
             ! neighbor 5 of particle's master
             proccoords(1) = nint( dims(1)*xcglob(idp)/lx - ax ) - 1
             proccoords(2) = nint( dims(2)*ycglob(idp)/ly - ay )
             call MPI_CART_RANK(comm_cart,proccoords,procrank,error)
             if (myid .eq. procrank) then
                if ( xcglob(idp) .lt. (leftbound+(radius+offset)) ) then
                   if(count_slve_loc.eq. 0 ) i = i+1
                   count_slve_loc = count_slve_loc + 1
                   ap(i)%mslv = -idp     ! myid is slave of particle abs(ap(p)%mslv)
                   ap(i)%nb(1) = 1      ! neighbor 1 of myid is particle's master
                endif
             endif
             ! neighbor 6 of particle's master
             proccoords(1) = nint( dims(1)*xcglob(idp)/lx - ax ) - 1
             proccoords(2) = nint( dims(2)*ycglob(idp)/ly - ay ) + 1
             call MPI_CART_RANK(comm_cart,proccoords,procrank,error)
             if (myid .eq. procrank) then
                dist = sqrt( (leftbound-xcglob(idp))**2. + (backbound-ycglob(idp))**2. )
                if ( abs(dist) .lt. (radius+offset) ) then
                   if(count_slve_loc.eq. 0 ) i = i+1
                   count_slve_loc = count_slve_loc + 1
                   ap(i)%mslv = -idp     ! myid is slave of particle abs(ap(p)%mslv)
                   ap(i)%nb(2) = 1      ! neighbor 2 of myid is particle's master
                endif
             endif
             ! neighbor 7 of particle's master
             proccoords(1) = nint( dims(1)*xcglob(idp)/lx - ax )
             proccoords(2) = nint( dims(2)*ycglob(idp)/ly - ay ) + 1
             call MPI_CART_RANK(comm_cart,proccoords,procrank,error)
             if (myid .eq. procrank) then
                if ( ycglob(idp) .gt. (backbound-(radius+offset)) ) then
                   if(count_slve_loc.eq. 0 ) i = i+1
                   count_slve_loc = count_slve_loc + 1
                   ap(i)%mslv = -idp     ! myid is slave of particle p=ap(p)%mslv
                   ap(i)%nb(3) = 1      ! neighbor 3 of myid is particle's master
                endif
             endif
             ! neighbor 8 of particle's master
             proccoords(1) = nint( dims(1)*xcglob(idp)/lx - ax ) + 1
             proccoords(2) = nint( dims(2)*ycglob(idp)/ly - ay ) + 1
             call MPI_CART_RANK(comm_cart,proccoords,procrank,error)
             if (myid .eq. procrank) then
                dist = sqrt( (rightbound-xcglob(idp))**2. + (backbound-ycglob(idp))**2. )
                if ( abs(dist) .lt. (radius+offset) ) then
                   if(count_slve_loc.eq. 0 ) i = i+1
                   count_slve_loc = count_slve_loc + 1
                   ap(i)%mslv = -idp     ! myid is slave of particle p=ap(p)%mslv
                   ap(i)%nb(4) = 1      ! neighbor 4 of myid is particle's master
                endif
             endif
          endif
       enddo
       pmax = i 
       !
       write(6,'(A7,I5,A8,I7,A18,I7,A11,A8,I5)') 'Thread ', myid, ' masters ', npmstr, ' and is slave for ', pmax-npmstr, ' particles. ', ' pmax = ', pmax
       !
       call MPI_ALLREDUCE(npmstr,count_mstr_all,1,MPI_INTEGER,MPI_SUM,comm_cart,error)
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
       else
          write(6,*) 'The particles were successfully initialized in thread ', myid, ' !'
       endif
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
                      call MPI_ISEND(ap(p)%x,11,MPI_REAL8,neighbor(nbsend), &
                           tag,comm_cart,arrayrequests((nrrequests-1) + 1),error)
                      ! send x,y,z,theta,phi,u,v,w,omx,omy,omz
                   endif
                endif
             endif
             if (ap(p)%mslv .lt. 0) then
                ! myid is slave of particle -ap(p)%mslv
                if (ap(p)%nb(nbrecv) .eq. 1) then
                   ! neighbor(nbrecv) is rank of master of particle -ap(p)%mslv
                   nrrequests = nrrequests + 1
                   call MPI_IRECV(ap(p)%x,11,MPI_REAL8,neighbor(nbrecv), &
                        tag,comm_cart,arrayrequests((nrrequests-1) + 1),error)
                   ! recv x,y,z,theta,phi,u,v,w,omx,omy,omz
                endif
             endif
          enddo ! do nb=
          call MPI_WAITALL(nrrequests,arrayrequests,arraystatuses,error)
       enddo
       !
       ! initial particle positions written to file (COBP)
       !
       !  write(rankpr,'(i4.4)') myid
       !  open(25,file=datadir//'mslv'//rankpr//'.txt')
       !  do p=1,pmax
       !    if (ap(p)%mslv .gt. 0) then
       !      counter = 0
       !      do i=1,8
       !        if (ap(p)%nb(i) .eq. 1) then
       !          counter = counter+1
       !          write(25,'(I4,A1,I5,A1,I5,A1,I2,A1,I2,2E16.8)') &
       !                myid,' ',p,' ',ap(p)%mslv,' ',i,' ',neighbor(i),ap(p)%x,ap(p)%y
       !          write(6,'(A29,I4,A1,I5,A1,I5,A1,I2,A1,I2,2E16.8)') 'rank,p,pms,nbr,ranknbr,x,y = ', &
       !                myid,' ',p,' ',ap(p)%mslv,' ',i,' ',neighbor(i),ap(p)%x,ap(p)%y
       !        endif
       !      enddo
       !      ! in case of no overlap with any neighbor
       !      if (counter .eq. 0) then
       !        write(25,'(I4,A1,I5,A1,I5,A1,I2,A1,I2,2E16.8)') &
       !                myid,' ',p,' ',ap(p)%mslv,' ',99,' ',99,ap(p)%x,ap(p)%y
       !        write(6,'(A29,I4,A1,I5,A1,I5,A1,I2,A1,I2,2E16.8)') 'rank,p,pms,nbr,ranknbr,x,y = ', &
       !                myid,' ',p,' ',ap(p)%mslv,' ',99,' ',99,ap(p)%x,ap(p)%y
       !      endif
       !    endif
       !    if (ap(p)%mslv .lt. 0) then
       !      do i=1,8
       !        if (ap(p)%nb(i) .eq. 1) then
       !          write(25,'(I4,A1,I5,A1,I5,A1,I2,A1,I2,2E16.8)') &
       !                myid,' ',p,' ',ap(p)%mslv,' ',i,' ',neighbor(i),ap(p)%x,ap(p)%y
       !          write(6,'(A29,I4,A1,I5,A1,I5,A1,I2,A1,I2,2E16.8)') 'rank,p,pms,nbr,ranknbr,x,y = ', &
       !                myid,' ',p,' ',ap(p)%mslv,' ',i,' ',neighbor(i),ap(p)%x,ap(p)%y
       !        endif
       !      enddo
       !    endif
       !  enddo
       !  close(25)
       !
       ! write tecplot output file with particle distribution
       !
       !  if (myid.eq.0) then
       !    open(25,file=datadir//'distr_particles.plt',form='formatted')
       !    write(25,*) 'TITLE="Tecplot Output"'
       !    write(25,*) 'VARIABLES= "X" "Y"'
       !    do p=1,pmax
       !      write(tempstr2,'(A,I2.2)') 'Particle ',p
       !      write(25,*) 'ZONE F=POINT T="', tempstr2, '" I=',250,' J=1'
       !      do i=1,250
       !        angle = 2.*pi*i/250.
       !        xp = glob(p)%x+radius*cos(angle)
       !        yp = glob(p)%y+radius*sin(angle)
       !          if (xp.lt.0) then
       !            xp = xp + lx
       !          elseif(xp.gt.lx) then
       !            xp = xp - lx
       !          endif
       !          if (yp.lt.0) then
       !            yp = yp + ly
       !          elseif(yp.gt.ly) then
       !            yp = yp - ly
       !          endif 
       !       write(25,*) xp/lx,yp/ly
       !      enddo
       !    enddo
       !    close(25) 
       !  endif
       ! (COBP)
       deallocate(glob)
    endif ! in.eq.0
    !
    if (in.eq.1) then
       !
       ! write particle related data directly in parallel with MPI-IO
       !
       allocate(glob(max(1,npmstr))) !allocate(glob(npmstr))
       npmstr_glob(:) = 0
       npmstr_glob(myid) = npmstr
       call MPI_ALLREDUCE(npmstr_glob(0),npmstr_glob_all(0),product(dims),MPI_INTEGER,MPI_SUM,comm_cart,error)
       mydisp = 0
       if(myid.ne.0) mydisp = sum(npmstr_glob_all(0:myid-1))
       i = 0
       do p=1,pmax
          if(ap(p)%mslv.gt.0) then
             idp = ap(p)%mslv
             i = i + 1
             glob(i)%x = ap(p)%x !WRITE THIS IN A MORE COMPACT WAY!
             glob(i)%y = ap(p)%y
             glob(i)%z = ap(p)%z
             glob(i)%theta = ap(p)%theta
             glob(i)%phi = ap(p)%phi
             glob(i)%u = ap(p)%u
             glob(i)%v = ap(p)%v
             glob(i)%w = ap(p)%w
             glob(i)%omx = ap(p)%omx
             glob(i)%omy = ap(p)%omy
             glob(i)%omz = ap(p)%omz
             glob(i)%omtheta = ap(p)%omtheta
             glob(i)%intu = ap(p)%intu
             glob(i)%intv = ap(p)%intv
             glob(i)%intw = ap(p)%intw
             glob(i)%intomx = ap(p)%intomx
             glob(i)%intomy = ap(p)%intomy
             glob(i)%intomz = ap(p)%intomz
             glob(i)%fxltot = ap(p)%fxltot
             glob(i)%fyltot = ap(p)%fyltot
             glob(i)%fzltot = ap(p)%fzltot
             ! restore variable below for restarting purpose             
             glob(i)%colfx = ap(p)%colfx
             glob(i)%colfy = ap(p)%colfy
             glob(i)%colfz = ap(p)%colfz
             glob(i)%coltx = ap(p)%coltx
             glob(i)%colty = ap(p)%colty
             glob(i)%coltz = ap(p)%coltz
             glob(i)%dx(1:nqmax) = ap(p)%dx(1:nqmax)
             glob(i)%dy(1:nqmax) = ap(p)%dy(1:nqmax)
             glob(i)%dz(1:nqmax) = ap(p)%dz(1:nqmax)
             glob(i)%dxt(1:nqmax) = ap(p)%dxt(1:nqmax)
             glob(i)%dyt(1:nqmax) = ap(p)%dyt(1:nqmax)
             glob(i)%dzt(1:nqmax) = ap(p)%dzt(1:nqmax)
             glob(i)%firstc(1:nqmax) = 1.*ap(p)%firstc(1:nqmax)
             glob(i)%qmax = 1.*ap(p)%qmax
             glob(i)%idp = 1.*idp
          endif
       enddo
       call MPI_FILE_OPEN(MPI_COMM_WORLD, 'allpartdata'//istepchar, &
            MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, error)
       filesize = 0_MPI_OFFSET_KIND
       call MPI_FILE_SET_SIZE(fh,filesize,error)  ! guarantee overwriting
       disp = mydisp*skip*lenr
       call MPI_FILE_SET_VIEW(fh, disp, MPI_REAL8,MPI_REAL8, 'native', & 
            MPI_INFO_NULL, error)
       call MPI_FILE_WRITE(fh,glob(1)%x,skip*npmstr,MPI_REAL8,MPI_STATUS_IGNORE,error)
       call MPI_FILE_CLOSE(fh,error)
       deallocate(glob)
    endif
    !
    return
  end subroutine loadpart
  !
end module mod_loadpart
