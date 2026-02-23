module mod_initmpi
  use decomp_2d
  use mod_param
  use mod_common_mpi
  implicit none
  private
  public initmpi, left, right, front, back, neighbor
  !integer :: left, right, front, back
  !integer :: neighbor(0:8)
contains
  !
  subroutine initmpi
    implicit none
    integer :: coordsleft(1:ndims) ,coordsright(1:ndims)
    integer :: coordsfront(1:ndims),coordsback(1:ndims)
    integer :: coordsneighbor(1:ndims)
    character(len=5) rankpr 
    !
    call MPI_INIT(error)
    !call MPI_COMM_RANK(MPI_COMM_WORLD,rankcw,error) -> not needed, I use 2decomp&fft routines*
    !
    periods(1) = .true.
    periods(2) = .true.
    periods(3) = .false.
    !

    !call MPI_COMM_SIZE(MPI_COMM_WORLD,Nproc,error)  ! -> not needed, I use 2decomp&fft routines*
    !write(*, *), dims, Nproc

    call decomp_2d_init(itot, jtot, ktot, dims(1), dims(2), periods)
    myid = nrank
    
    coords(1) = (zstart(1)-1)*dims(1)/itot ! iglob = i + zstart(1) - 1 (SS)
    coords(2) = (zstart(2)-1)*dims(2)/jtot ! jglob = j + zstart(2) - 1 (SS)

   !  !
   !  ! Debug - output boundleftmyid value (SS_debug)
   !  if (myid == 10) then
   !     write(6,'(A,I0,A,E16.8)') 'Process ', myid, ': boundleftmyid = ', boundleftmyid
   !     write(6,'(A,I0,A,E16.8)') 'Process ', myid, ': boundfrontmyid = ', boundfrontmyid
   !     write(6,'(A,I0,A,E16.8)') 'Process ', myid, ': coords(1) = ', coords(1)
   !     write(6,'(A,I0,A,E16.8)') 'Process ', myid, ': coords(2) = ', coords(2)
   !     write(6,'(A,I0,A,E16.8)') 'Process ', myid, ': zstart(1) = ', zstart(1)
   !     write(6,'(A,I0,A,E16.8)') 'Process ', myid, ': zstart(2) = ', zstart(2)
   !  endif
    if ( jtot .lt. itot) then
       if (myid.eq.0) then
          write(6,*) 'jtot smaller than itot: this is computationally less efficient!'
          write(6,*) 'Interchange the coordinate directions.'
          write(6,*) 'Program aborted...'
       endif
       call mpi_finalize(error)
       stop
    endif
    if ( mod(itot,dims(1)).ne.0 ) then
       if (myid.eq.0) then
          write(6,*) 'itot not divisable by the number of threads in its direction'
          write(6,*) 'Change the value of itot or dims(1) in "param.f90".'
          write(6,*) 'Program aborted...'
       endif
       call mpi_finalize(error)
       stop
    endif
    if ( mod(jtot,dims(2)).ne.0 ) then
       if (myid.eq.0) then
          write(6,*) 'jtot not divisable by the number of threads in its direction'
          write(6,*) 'Change the value of jtot or dims(2) in "param.f90".'
          write(6,*) 'Program aborted...'
       endif
       call mpi_finalize(error)
       stop
    endif

    
    !
    !neighbor(6)/leftback        neighbor(7)/back     neighbor(8)/rightback
    !
    !neighbor(5)/left            neighbor(0)/myid     neighbor(1)/right
    !
    !neighbor(4)/leftfront       neighbor(3)/front    neighbor(2)/rightfront
    !
    !create virtual topology: Cartesian grid *not needed, I use 2decomp&fft routines*
    ! periods(1) = .true. !periodic b.c.'s in x-dir.
    ! periods(2) = .true. !periodic b.c.'s in y-dir.
    ! reorder    = .true. !reorder ranks 
    ! call MPI_CART_CREATE(MPI_COMM_WORLD,ndims,dims,periods(1:2), &
    !                      reorder,comm_cart,error)
    !myid is rank within group of communicator comm_cart
    ! call MPI_COMM_RANK  (comm_cart,myid,error)
    !
    comm_cart = DECOMP_2D_COMM_CART_Z
    !
    call MPI_CART_SHIFT(comm_cart,0,1,left,right,error)
    call MPI_CART_SHIFT(comm_cart,1,1,front,back,error)
    call MPI_CART_COORDS(comm_cart,right,ndims,coordsright,error)
    call MPI_CART_COORDS(comm_cart,front,ndims,coordsfront,error)
    coordsneighbor(1) = coordsright(1)
    coordsneighbor(2) = coordsfront(2)
    call MPI_CART_RANK(comm_cart,coordsneighbor,rightfront,error)
    call MPI_CART_COORDS(comm_cart,back ,ndims,coordsback ,error)
    coordsneighbor(1) = coordsright(1)
    coordsneighbor(2) = coordsback(2)
    call MPI_CART_RANK(comm_cart,coordsneighbor,rightback,error)
    call MPI_CART_COORDS(comm_cart,left,ndims,coordsleft,error)
    call MPI_CART_COORDS(comm_cart,front,ndims,coordsfront,error)
    coordsneighbor(1) = coordsleft(1)
    coordsneighbor(2) = coordsfront(2)
    call MPI_CART_RANK(comm_cart,coordsneighbor,leftfront,error)
    call MPI_CART_COORDS(comm_cart,back ,ndims,coordsback ,error)
    coordsneighbor(1) = coordsleft(1)
    coordsneighbor(2) = coordsback(2)
    call MPI_CART_RANK(comm_cart,coordsneighbor,leftback,error)
    !
    neighbor(0) = myid
    neighbor(1) = right
    neighbor(2) = rightfront
    neighbor(3) = front
    neighbor(4) = leftfront
    neighbor(5) = left
    neighbor(6) = leftback
    neighbor(7) = back
    neighbor(8) = rightback
    !
    !write(6,*) 'Neighbors of process ',myid,' : '
    !write(6,*)  leftback,back,rightback
    !write(6,*)  left,myid,right
    !write(6,*)  leftfront,front,rightfront
    !write(rankpr,'(i5.5)') myid
    !open(25,file=datadir//'topology'//rankpr//'.txt')
    !  write(25,*) 'Neighbors of process ',myid,' : '
    !  write(25,*)  leftback,back,rightback
    !  write(25,*)  left,myid,right
    !  write(25,*)  leftfront,front,rightfront
    !close(25)
    !
    ! Definitions of datatypes for velocity and pressure b.c.'s
    ! Note: array(i,j,k) is basically a 1-dim array; 
    !       k is the outer and i is the inner loop counter =>
    !         * for fixed i, (j1+1)*(k1+1) blocks of 1 element,
    !           with (i1+1) elements between start and end 
    !         * for fixed j, (k1+1) blocks of (i1+1) elements,
    !           with (i1+1)*(j1+1) elements between start and end
    !
    !
    call MPI_TYPE_VECTOR((j1+1)*(k1+1),1,(i1+1),MPI_REAL8,xhalo,error) 
    call MPI_TYPE_COMMIT(xhalo,error)
    call MPI_TYPE_VECTOR((k1+1),(i1+1),(i1+1)*(j1+1),MPI_REAL8,yhalo,error) 
    call MPI_TYPE_COMMIT(yhalo,error)
    !
    call MPI_TYPE_VECTOR((j1+1+2)*(k1+1+2),2,(i1+1+2),MPI_REAL8,xhalo_ibm,error)
    call MPI_TYPE_COMMIT(xhalo_ibm,error)
    call MPI_TYPE_VECTOR((k1+1+2),2*(i1+1+2),(i1+1+2)*(j1+1+2),MPI_REAL8,yhalo_ibm,error)
    call MPI_TYPE_COMMIT(yhalo_ibm,error)
    !
    return
  end subroutine initmpi
  !
end module mod_initmpi
