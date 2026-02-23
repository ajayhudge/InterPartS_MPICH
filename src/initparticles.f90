module mod_initparticles
  use decomp_2d
  use mod_param
  use mod_common
  use mod_common_mpi
  implicit none
  private
  public initparticles
contains
  subroutine initparticles
    implicit none
    !integer,parameter extrapoints = nint(2*(radius+0.1)*dxi)! for visualisation purposes
    real, dimension(np) :: xcglob,ycglob,zcglob,thetacglob,phicglob
    integer :: i,j,k,p,pp,rk
    integer :: proccoords(1:ndims),procrank
    integer, dimension(2) :: sbuf,rbuf
    real :: leftbound,rightbound,frontbound,backbound
    real :: dist,distx,disty,distz,distzw,angle
    real :: xp,yp
    real :: ax
    real :: ay
    real :: rn
    integer :: counter,crys
    integer :: idp
    character(len=5) rankpr 
    character(len=7) :: tempstr
    character(len=11) :: tempstr2
    integer :: count_mstr,count_mstr_all,count_slve_loc
    integer(kind=8), allocatable, dimension(:) :: iseed
    !
    ! position of spheres: global initialisation by root (myid = 0)
    !
    allocate(iseed(1))
    iseed(1) = 16112006
    crys=0
    if (myid .eq. 0) then
       call random_seed(put = iseed)
       write(6,*) 'part as crystal 1, random 0: ',crys
       open(23,file="position_spheres.txt")
       if(crys.eq.1) then
          !    p=0 ! (COBP) in the next 20 lines
          !    ! pseudo-crystal
          !    do k=1,4
          !      do j=1,16
          !        do i=1,8
          !          p=p+1
          !          thetacglob(p) = 0.
          !          phicglob(p)   = 0.
          !          call random_number(rn)
          !          xcglob(p)     = (lx/4.)*((1.*i)-0.5) + 0.25*(rn-0.5)*radius
          !          call random_number(rn)
          !          ycglob(p)     = (ly/8.)*((1.*j)-0.5) + 0.25*(rn-0.5)*radius
          !          call random_number(rn)
          !          zcglob(p)     = (lz/2. )*((1.*k)-0.5) + (rn-0.5)*radius
          !          write(6,*) 'Location of sphere #',p
          !          write(6,*) 'x,y,z = ',xcglob(p),ycglob(p),zcglob(p)
          !          write(23,'(I3,5E32.16)') p,thetacglob(p),phicglob(p), &
          !                                  xcglob(p),ycglob(p),zcglob(p)
          !        enddo
          !      enddo
          !    enddo
       else
          ! pseudo-random
          counter = 0
          do p=1,NP
             if(NP.eq.1) then
                thetacglob(p) = 0.
                phicglob(p)   = 0.
                xcglob(p)     = lx*0.5
                ycglob(p)     = ly*0.5
                zcglob(p)     = lz*0.9
             else      
111             continue
                thetacglob(p) = 0.
                phicglob(p)   = 0.
                call random_number(rn)
                xcglob(p)     = lx*rn
                call random_number(rn)
                ycglob(p)     = ly*rn
                call random_number(rn)
                zcglob(p)     = lz*rn
                distzw = min(abs(lz-zcglob(p)),abs(zcglob(p)))
                if(distzw.lt.1.05*radius) goto 111
                do pp=1,p
                   if (pp.eq.p) goto 444 ! could be changed by changing the loop limits!
                   distz = abs(zcglob(p)- zcglob(pp))
                   do j=-1,1
                      disty = abs(ycglob(p)- (ycglob(pp)+j*ly))
                      if(disty.gt.2.05*radius) goto 222
                      do i=-1,1
                         distx = abs(xcglob(p)- (xcglob(pp)+i*lx))
                         if(distx.gt.2.05*radius) goto 333
                         if(distx.gt.2.05*radius.or. &
                              disty.gt.2.05*radius.or. &
                              distz.gt.2.05*radius.or.p.eq.pp) then
                            ! good particle
                         else
                            dist = distx**2+disty**2.+distz**2.
                            if((dist.lt.(4.2*radius**2.))) then
                               !write(*,*)'RANDOM DEVIATION'
                               !write(*,*)dist,distw
                               counter=counter+1
                               goto 111
                            endif
                         endif
333                      continue
                      enddo
222                   continue
                   enddo
444                continue
                enddo
             endif
             write(6,*) 'Location of sphere #',p
             write(6,*) 'x,y,z = ',xcglob(p),ycglob(p),zcglob(p)
             write(23,'(I6,5E32.16)') p,thetacglob(p),phicglob(p), &
                  xcglob(p),ycglob(p),zcglob(p)
          enddo
          write(*,*)'RANDOM DEVIATIONS: ',counter
          p=p-1
          !    p = np
       endif
       close(23)
       if ( (p.ne.np) ) then
          print*,counter,np
          write(6,*) 'Fatal error in initialisation of particle positions!'
          write(6,*) 'Program aborted...'
          call mpi_finalize(error)
          stop
       endif
       do rk=1,Nproc-1
          call MPI_SSEND(xcglob    ,np,MPI_REAL8,rk,rk+0*(Nproc-1),comm_cart,error)
          call MPI_SSEND(ycglob    ,np,MPI_REAL8,rk,rk+1*(Nproc-1),comm_cart,error)
          call MPI_SSEND(zcglob    ,np,MPI_REAL8,rk,rk+2*(Nproc-1),comm_cart,error)
          call MPI_SSEND(thetacglob,np,MPI_REAL8,rk,rk+3*(Nproc-1),comm_cart,error)
          call MPI_SSEND(phicglob  ,np,MPI_REAL8,rk,rk+4*(Nproc-1),comm_cart,error)
       enddo
    else ! if myid is not 0:
       call MPI_RECV(xcglob    ,np,MPI_REAL8,0,myid+0*(Nproc-1),comm_cart,status,error)
       call MPI_RECV(ycglob    ,np,MPI_REAL8,0,myid+1*(Nproc-1),comm_cart,status,error)
       call MPI_RECV(zcglob    ,np,MPI_REAL8,0,myid+2*(Nproc-1),comm_cart,status,error)
       call MPI_RECV(thetacglob,np,MPI_REAL8,0,myid+3*(Nproc-1),comm_cart,status,error)
       call MPI_RECV(phicglob  ,np,MPI_REAL8,0,myid+4*(Nproc-1),comm_cart,status,error)
    endif
    !
    ! Determine master and slave processes for each particle.
    !
    ! initialisation
    !
    ap(1:npmax)%x = 0.
    ap(1:npmax)%y = 0.
    ap(1:npmax)%z = 0.
    ap(1:npmax)%theta = 0.
    ap(1:npmax)%phi = 0.
    ap(1:npmax)%mslv = 0.
    forall(i=1:npmax) ap(i)%nb(1:8) = 0
    !
    count_mstr = 0
    i = 0
    ax = 0.5
    ay = 0.5
    !
    pmax = 0
    do p=1,np
       if (xcglob(p).lt.0..or.xcglob(p).gt.lx .or. &
            ycglob(p).lt.0..or.ycglob(p).gt.ly .or. &
            zcglob(p).lt.radius.or.zcglob(p).gt.lz-radius) then
          if (myid.eq.0) then
             write(6,*) 'Fatal error in initialisation of particle positions - '
             write(6,*) 'particle outside the domain!'
             write(6,*) 'Program aborted...'
          endif
          call mpi_finalize(error)
          stop
       endif
       if (xcglob(p).eq.lx) ax = 0.51
       if (xcglob(p).eq.0) ax = 0.49
       if (ycglob(p).eq.ly) ay = 0.51
       if (ycglob(p).eq.0) ay = 0.49
       !  
       proccoords(1) = nint(xcglob(p)*dims(1)/lx - ax)
       proccoords(2) = nint(ycglob(p)*dims(2)/ly - ay)
       leftbound     = (proccoords(1)  )*lx/(1.*dims(1)) ! left  boundary of particle's master
       rightbound    = (proccoords(1)+1)*lx/(1.*dims(1)) ! right boundary of particle's master
       frontbound    = (proccoords(2)  )*ly/(1.*dims(2)) ! front boundary of particle's master
       backbound     = (proccoords(2)+1)*ly/(1.*dims(2)) ! back  boundary of particle's master
       call MPI_CART_RANK(comm_cart,proccoords,procrank,error)
       if (myid.eq.procrank) then
          i = i + 1
          count_mstr = count_mstr + 1
          ap(i)%x = xcglob(p)
          ap(i)%y = ycglob(p)
          ap(i)%z = zcglob(p)
          ap(i)%theta = thetacglob(p)
          ap(i)%phi = phicglob(p)
          ap(i)%mslv = p
          !neighbor 1
          if ( ap(i)%x .gt. (rightbound-(radius+offset))) then 
             ap(i)%nb(1) = 1 !neighbor 1 is slave of particle ap(i)%mslv 
          endif
          !neighbor 2
          dist = sqrt( (rightbound-ap(i)%x)**2. + (frontbound-ap(i)%y)**2. )
          if ( abs(dist) .lt. (radius+offset) ) then
             ap(i)%nb(2) = 1 !neighbor 2 is slave of particle ap(i)%mslv
          endif
          !neighbor 3
          if ( ap(i)%y .lt. (frontbound+(radius+offset))) then
             ap(i)%nb(3) = 1 !neighbor 3 is slave of particle ap(i)%mslv
          endif
          !neighbor 4
          dist = sqrt( (leftbound-ap(i)%x)**2. + (frontbound-ap(i)%y)**2. ) 
          if ( abs(dist) .lt. (radius+offset)) then
             ap(i)%nb(4) = 1 !neighbor 4 is slave of particle ap(i)%mslv
          endif
          !neighbor 5
          if ( ap(i)%x .lt. (leftbound+(radius+offset)) ) then
             ap(i)%nb(5) = 1 !neighbor 5 is slave of particle ap(i)%mslv
          endif
          !neighbor 6
          dist = sqrt( (leftbound-ap(i)%x)**2. + (backbound-ap(i)%y)**2. )
          if ( abs(dist) .lt. (radius+offset) ) then
             ap(i)%nb(6) = 1 !neighbor 6 is slave of particle ap(i)%mslv
          endif
          !neighbor 7
          if ( ap(i)%y .gt. (backbound-(radius+offset)) ) then
             ap(i)%nb(7) = 1 !neighbor 7 is slave of particle ap(i)%mslv
          endif
          !neighbor 8
          dist = sqrt( (rightbound-ap(i)%x)**2. + (backbound-ap(i)%y)**2. )
          if ( abs(dist) .lt. (radius+offset) ) then
             ap(i)%nb(8) = 1 !neighbor 8 is slave of particle ap(i)%mslv
          endif
       else
          count_slve_loc = 0
          !neighbor 1 of particle's master
          proccoords(1) = nint( dims(1)*xcglob(p)/lx - ax ) + 1
          proccoords(2) = nint( dims(2)*ycglob(p)/ly - ay ) 
          call MPI_CART_RANK(comm_cart,proccoords,procrank,error)
          if (myid .eq. procrank) then
             if ( xcglob(p) .gt. (rightbound-(radius+offset))) then
                if(count_slve_loc.eq.0) i = i+1
                ap(i)%mslv = -p     !myid is slave of particle abs(ap(i)%mslv)
                ap(i)%nb(5) = 1      !neighbor 5 of myid is particle's master
                count_slve_loc = count_slve_loc + 1
                ap(i)%x = xcglob(p)
                ap(i)%y = ycglob(p)
                ap(i)%z = zcglob(p)
                ap(i)%theta = thetacglob(p)
                ap(i)%phi = phicglob(p)
             endif
          endif
          !neighbor 2 of particle's master
          proccoords(1) = nint( dims(1)*xcglob(p)/lx - ax ) + 1
          proccoords(2) = nint( dims(2)*ycglob(p)/ly - ay ) - 1
          call MPI_CART_RANK(comm_cart,proccoords,procrank,error)
          if (myid .eq. procrank) then
             dist = sqrt( (rightbound-xcglob(p))**2. + (frontbound-ycglob(p))**2. ) 
             if ( abs(dist) .lt. (radius+offset) ) then
                if(count_slve_loc.eq.0) i = i+1
                ap(i)%mslv = -p     !myid is slave of particle abs(ap(i)%mslv)
                ap(i)%nb(6) = 1      !neighbor 6 of myid is particle's master
                count_slve_loc = count_slve_loc + 1
                ap(i)%x = xcglob(p)
                ap(i)%y = ycglob(p)
                ap(i)%z = zcglob(p)
                ap(i)%theta = thetacglob(p)
                ap(i)%phi = phicglob(p)
             endif
          endif
          !neighbor 3 of particle's master
          proccoords(1) = nint( dims(1)*xcglob(p)/lx - ax )
          proccoords(2) = nint( dims(2)*ycglob(p)/ly - ay ) - 1
          call MPI_CART_RANK(comm_cart,proccoords,procrank,error)
          if (myid .eq. procrank) then
             if ( ycglob(p) .lt. (frontbound+(radius+offset)) ) then
                if(count_slve_loc.eq.0) i = i+1
                ap(i)%mslv = -p     !myid is slave of particle abs(ap(i)%mslv)
                ap(i)%nb(7) = 1      !neighbor 7 of myid is particle's master
                count_slve_loc = count_slve_loc + 1
                ap(i)%x = xcglob(p)
                ap(i)%y = ycglob(p)
                ap(i)%z = zcglob(p)
                ap(i)%theta = thetacglob(p)
                ap(i)%phi = phicglob(p)
             endif
          endif
          !neighbor 4 of particle's master
          proccoords(1) = nint( dims(1)*xcglob(p)/lx - ax ) - 1
          proccoords(2) = nint( dims(2)*ycglob(p)/ly - ay ) - 1
          call MPI_CART_RANK(comm_cart,proccoords,procrank,error)
          if (myid .eq. procrank) then
             dist = sqrt( (leftbound-xcglob(p))**2. + (frontbound-ycglob(p))**2. )
             if ( abs(dist) .lt. (radius+offset) ) then
                if(count_slve_loc.eq.0) i = i+1
                ap(i)%mslv = -p     !myid is slave of particle abs(ap(i)%mslv)
                ap(i)%nb(8) = 1      !neighbor 8 of myid is particle's master
                count_slve_loc = count_slve_loc + 1
                ap(i)%x = xcglob(p)
                ap(i)%y = ycglob(p)
                ap(i)%z = zcglob(p)
                ap(i)%theta = thetacglob(p)
                ap(i)%phi = phicglob(p)
             endif
          endif
          !neighbor 5 of particle's master
          proccoords(1) = nint( dims(1)*xcglob(p)/lx - ax ) - 1
          proccoords(2) = nint( dims(2)*ycglob(p)/ly - ay )
          call MPI_CART_RANK(comm_cart,proccoords,procrank,error)
          if (myid .eq. procrank) then
             if ( xcglob(p) .lt. (leftbound+(radius+offset)) ) then
                if(count_slve_loc.eq.0) i = i+1
                ap(i)%mslv = -p     !myid is slave of particle abs(ap(i)%mslv)
                ap(i)%nb(1) = 1      !neighbor 1 of myid is particle's master
                count_slve_loc = count_slve_loc + 1
                ap(i)%x = xcglob(p)
                ap(i)%y = ycglob(p)
                ap(i)%z = zcglob(p)
                ap(i)%theta = thetacglob(p)
                ap(i)%phi = phicglob(p)
             endif
          endif
          !neighbor 6 of particle's master
          proccoords(1) = nint( dims(1)*xcglob(p)/lx - ax ) - 1
          proccoords(2) = nint( dims(2)*ycglob(p)/ly - ay ) + 1
          call MPI_CART_RANK(comm_cart,proccoords,procrank,error)
          if (myid .eq. procrank) then
             dist = sqrt( (leftbound-xcglob(p))**2. + (backbound-ycglob(p))**2. )
             if ( abs(dist) .lt. (radius+offset) ) then
                if(count_slve_loc.eq.0) i = i+1
                ap(i)%mslv = -p     !myid is slave of particle abs(ap(i)%mslv)
                ap(i)%nb(2) = 1      !neighbor 2 of myid is particle's master
                count_slve_loc = count_slve_loc + 1
                ap(i)%x = xcglob(p)
                ap(i)%y = ycglob(p)
                ap(i)%z = zcglob(p)
                ap(i)%theta = thetacglob(p)
                ap(i)%phi = phicglob(p)
             endif
          endif
          !neighbor 7 of particle's master
          proccoords(1) = nint( dims(1)*xcglob(p)/lx - ax )
          proccoords(2) = nint( dims(2)*ycglob(p)/ly - ay ) + 1
          call MPI_CART_RANK(comm_cart,proccoords,procrank,error)
          if (myid .eq. procrank) then
             if ( ycglob(p) .gt. (backbound-(radius+offset)) ) then
                if(count_slve_loc.eq.0) i = i+1
                ap(i)%mslv = -p     !myid is slave of particle p=ap(i)%mslv
                ap(i)%nb(3) = 1      !neighbor 3 of myid is particle's master
                count_slve_loc = count_slve_loc + 1
                ap(i)%x = xcglob(p)
                ap(i)%y = ycglob(p)
                ap(i)%z = zcglob(p)
                ap(i)%theta = thetacglob(p)
                ap(i)%phi = phicglob(p)
             endif
          endif
          !neighbor 8 of particle's master
          proccoords(1) = nint( dims(1)*xcglob(p)/lx - ax ) + 1
          proccoords(2) = nint( dims(2)*ycglob(p)/ly - ay ) + 1
          call MPI_CART_RANK(comm_cart,proccoords,procrank,error)
          if (myid .eq. procrank) then
             dist = sqrt( (rightbound-xcglob(p))**2. + (backbound-ycglob(p))**2. )
             if ( abs(dist) .lt. (radius+offset) ) then
                if(count_slve_loc.eq.0) i = i+1
                ap(i)%mslv = -p     !myid is slave of particle p=ap(i)%mslv
                ap(i)%nb(4) = 1      !neighbor 4 of myid is particle's master
                count_slve_loc = count_slve_loc + 1
                ap(i)%x = xcglob(p)
                ap(i)%y = ycglob(p)
                ap(i)%z = zcglob(p)
                ap(i)%theta = thetacglob(p)
                ap(i)%phi = phicglob(p)
             endif
          endif
       endif
    enddo
    !
    ! maximum number of particles in a thread is equal to the number of particles
    ! 'mastered' and 'slaved' by it
    !
    pmax = i 
    npmstr = count_mstr
    write(6,'(A7,I5,A8,I5,A18,I5,A11,A8,I5)') 'Thread ', myid, ' masters ', count_mstr, ' and is slave for ', pmax-npmstr, ' particles. ', ' pmax = ', pmax
    !
    call MPI_ALLREDUCE(count_mstr,count_mstr_all,1,MPI_INTEGER,MPI_SUM,comm_cart,error)
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
    !
    ! initial particle positions written to file
    !
    !write(rankpr,'(i5.5)') myid
    !open(25,file=datadir//'mslv'//rankpr//'.txt')
    !do p=1,pmax
    !  idp = abs(ap(p)%mslv)
    !  if (ap(p)%mslv .gt. 0) then
    !    counter = 0
    !    do i=1,8
    !      if (ap(p)%nb(i) .eq. 1) then
    !        counter = counter+1
    !        write(25,'(I4,A1,I5,A1,I5,A1,I2,A1,I5,2E16.8)') &
    !              myid,' ',idp,' ',ap(p)%mslv,' ',i,' ',neighbor(i),ap(p)%x,ap(p)%y
    !        write(6,'(A29,I4,A1,I5,A1,I5,A1,I2,A1,I5,2E16.8)') 'rank,p,pms,nbr,ranknbr,x,y = ', &
    !              myid,' ',idp,' ',ap(p)%mslv,' ',i,' ',neighbor(i),ap(p)%x,ap(p)%y
    !      endif
    !    enddo
    !    ! in case of no overlap with any neighbor
    !    if (counter .eq. 0) then
    !      write(25,'(I4,A1,I5,A1,I5,A1,I2,A1,I5,2E16.8)') &
    !              myid,' ',idp,' ',ap(p)%mslv,' ',99,' ',99,ap(p)%x,ap(p)%y
    !      write(6,'(A29,I4,A1,I5,A1,I5,A1,I2,A1,I5,2E16.8)') 'rank,p,pms,nbr,ranknbr,x,y = ', &
    !              myid,' ',idp,' ',ap(p)%mslv,' ',99,' ',99,ap(p)%x,ap(p)%y
    !    endif
    !  endif
    !  if (ap(p)%mslv .lt. 0) then
    !    do i=1,8
    !      if (ap(p)%nb(i) .eq. 1) then
    !        write(25,'(I4,A1,I5,A1,I5,A1,I2,A1,I5,2E16.8)') &
    !              myid,' ',idp,' ',ap(p)%mslv,' ',i,' ',neighbor(i),ap(p)%x,ap(p)%y
    !        write(6,'(A29,I4,A1,I5,A1,I5,A1,I2,A1,I5,2E16.8)') 'rank,p,pms,nbr,ranknbr,x,y = ', &
    !              myid,' ',idp,' ',ap(p)%mslv,' ',i,' ',neighbor(i),ap(p)%x,ap(p)%y
    !      endif
    !    enddo
    !  endif
    !enddo
    !close(25)
    !!
    !! write topology of the domain decomposition in a tecplot file
    !!
    !if (myid==0) then
    !  open(25,file=datadir//'decomp.plt',form='formatted')
    !  write(25,*) 'TITLE="Tecplot Output"'
    !  write(25,*) 'VARIABLES= "X" "Y" "VAR"'
    !  write(25,*) 'ZONE F=POINT T="Rank 00" I=',imax,' J=',jmax
    !  do i=1,imax
    !    do j=1,jmax
    !      write(25,*) (i+zstart(1)-1)*dx/lx,(j+zstart(2)-1)*dy/ly,myid
    !    enddo
    !  enddo
    !  do rk=1,nproc-1
    !    CALL MPI_RECV(rbuf,2,MPI_INTEGER,rk,rk,MPI_COMM_WORLD,status,error)
    !    write(tempstr,'(A,I2.2)') 'Rank ',rk
    !  write(25,*) 'TITLE="Tecplot Output"'
    !  write(25,*) 'VARIABLES= "X" "Y" "VAR"'
    !  write(25,*) 'ZONE F=POINT T="', tempstr, '" I=',imax,' J=',jmax
    !    do i=1,imax
    !      do j=1,jmax
    !        write(25,*) (i+rbuf(1)-1)*dx/lx,(j+rbuf(2)-1)*dy/ly,rk
    !      enddo
    !    enddo
    !  enddo
    !  close(25)
    !else
    !  sbuf(1) = zstart(1)
    !  sbuf(2) = zstart(2)
    !  CALL MPI_SEND(sbuf,2,MPI_INTEGER,0,myid,MPI_COMM_WORLD,error)
    !endif
    !
    ! write initial particle distribution in a tecplot file
    !
    !if (myid.eq.0) then
    !  open(25,file=datadir//'distr_particles.plt',form='formatted')
    !  write(25,*) 'TITLE="Tecplot Output"'
    !  write(25,*) 'VARIABLES= "X" "Y"'
    !  do p=1,np
    !    write(tempstr2,'(A,I2.2)') 'Particle ',p
    !    write(25,*) 'ZONE F=POINT T="', tempstr2, '" I=',250,' J=1'
    !    do i=1,250
    !      angle = 2.*pi*i/250.
    !      xp = xcglob(p)+radius*cos(angle)
    !      yp = ycglob(p)+radius*sin(angle)
    !        if (xp.lt.0) then
    !          xp = xp + lx
    !        elseif(xp.gt.lx) then
    !          xp = xp - lx
    !        endif
    !        if (yp.lt.0) then
    !          yp = yp + ly
    !        elseif(yp.gt.ly) then
    !          yp = yp - ly
    !        endif 
    !      write(25,*) xp/lx,yp/ly
    !    enddo
    !  enddo
    !  close(25) 
    !endif
    !
    ap(1:npmax)%u = 0.
    ap(1:npmax)%v = 0.
    ap(1:npmax)%w = 0.
    ap(1:npmax)%omx = 0.
    ap(1:npmax)%omy = 0.
    ap(1:npmax)%omz = 0.
    ap(1:npmax)%omtheta = 0.
    ap(1:npmax)%vol = volp
    ap(1:npmax)%mominert = mominert
    ap(1:npmax)%ratiorho = ratiorho
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
    end forall
    !
    return
  end subroutine initparticles
  !
end module mod_initparticles
