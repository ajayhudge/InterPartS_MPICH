module mod_bound
  use mod_common
  use decomp_2d
  use mod_param
  use mod_common_mpi
  implicit none
  private
  public bounduvw, boundp, updthalos,updthalos_ibm, boundscal
contains
  !
  subroutine bounduvw(u,v,w)
    use mod_param
    implicit none
    integer :: i,j,k
    real, dimension(0:,0:,0:) :: u,v,w

    select case(bc_uvw_top_type)
    case('BC_NOSLIP')
          !$omp parallel default(none) shared(u,v,w) private(i,j)
          !$omp do
          do j=0,j1
             do i=0,i1
                u(i,j,k1)   = -u(i,j,kmax)    ! no-slip
                v(i,j,k1)   = -v(i,j,kmax)    ! no-slip
                w(i,j,kmax) = 0.             ! no-penetration
             enddo
          enddo
          !$omp end parallel
    case('BC_FREESLIP')
          !$omp parallel default(none) shared(u,v,w) private(i,j)
          !$omp do
          do j=0,j1
             do i=0,i1
                u(i,j,k1)   = u(i,j,kmax)    ! free-slip
                v(i,j,k1)   = v(i,j,kmax)    ! free-slip
                w(i,j,kmax) = 0.             ! no-penetration
             enddo
          enddo
          !$omp end parallel
    end select

    select case(bc_uvw_bot_type)
    case('BC_NOSLIP')
          !$omp parallel default(none) shared(u,v,w) private(i,j)
          !$omp do
          do j=0,j1
             do i=0,i1
                u(i,j,0)   = -u(i,j,1)    ! no-slip
                v(i,j,0)   = -v(i,j,1)    ! no-slip
                w(i,j,0)   = 0.             ! no-penetration
             enddo
          enddo
          !$omp end parallel
    case('BC_FREESLIP')
          !$omp parallel default(none) shared(u,v,w) private(i,j)
          !$omp do
          do j=0,j1
             do i=0,i1
                u(i,j,0)   = u(i,j,1)    ! free-slip
                v(i,j,0)   = v(i,j,1)    ! free-slip
                w(i,j,0)   = 0.          ! no-penetration
             enddo
          enddo
          !$omp end parallel
    end select

    !
    ! communicate data in x direction (periodic b.c.'s incorporated)
    call updthalos(u,1)
    call updthalos(v,1)
    call updthalos(w,1)
    !
    ! communicate data in y direction (periodic b.c.'s incorporated)
    !
    call updthalos(u,2)
    call updthalos(v,2)
    call updthalos(w,2)
    !
    return
  end subroutine bounduvw
  !
  subroutine boundp(p)
    use mod_param
    implicit none
    integer :: i,j,k
    real, dimension(0:,0:,0:) :: p
    !
    !$omp parallel default(none) shared(p) private(i,j)
    !$omp do
    do j=0,j1
       do i=0,i1
          p(i,j,0) = p(i,j,1)     ! Neumann (consistent with no/free-slip)
          p(i,j,k1) = p(i,j,kmax) ! Neumann (consistent with no/free-slip)
       enddo
    enddo
    !$omp end parallel

    !
    ! communicate data in x direction (periodic b.c.'s incorporated)
    call updthalos(p,1)
    !
    !  communicate data in y direction (periodic b.c.'s incorporated)
    !
    call updthalos(p,2)
    !
    return
  end subroutine boundp
  !
  subroutine boundscal(c)
    use mod_param
    implicit none
    integer :: i,j,k
    real, dimension(0:,0:,0:) :: c

    select case(bc_c_top_type)
    case('BC_DIRICHLET') 
    !$omp parallel default(none) shared(c) private(i,j)
    !$omp do
    do j=0,j1
       do i=0,i1
          c(i,j,k1) = 2*bc_c_top_val - c(i,j,kmax)  ! Dirichlet
       enddo
    enddo
    !$omp end parallel

    case('BC_NEUMANN') 
    !$omp parallel default(none) shared(c) private(i,j)
    !$omp do
    do j=0,j1
       do i=0,i1
          c(i,j,k1) = c(i,j,kmax) ! Neumann
       enddo
    enddo
    !$omp end parallel 
    end select     

    select case(bc_c_bot_type)
    case('BC_DIRICHLET') 
    !$omp parallel default(none) shared(c) private(i,j)
    !$omp do
    do j=0,j1
       do i=0,i1
          c(i,j,0) = 2*bc_c_bot_val - c(i,j,1)  ! Dirichlet
       enddo
    enddo
    !$omp end parallel
    
    case('BC_NEUMANN') 
    !$omp parallel default(none) shared(c) private(i,j)
    !$omp do
    do j=0,j1
       do i=0,i1
          c(i,j,0) = c(i,j,1)     ! Neumann
       enddo
    enddo
    !$omp end parallel 
    end select       

    !
    ! communicate data in x direction (periodic b.c.'s incorporated)
    !
    call updthalos(c,1)
    !
    ! communicate data in y direction (periodic b.c.'s incorporated)
    !
    call updthalos(c,2)
    !

    return

  end subroutine boundscal
  !
  subroutine updthalos(var,dir)
    use mpi
    use mod_param
    use mod_common_mpi
    implicit none
    real, dimension(0:,0:,0:), intent(inout) :: var
    integer, intent(in) :: dir
    integer :: requests(4), statuses(MPI_STATUS_SIZE,4)
    !
    !  This subroutine updates the halos that store info
    !  from the neighboring computational sub-domain
    !
    select case(dir)
    case(1) ! x direction
       call MPI_SENDRECV(var(1,0,0),1,xhalo,left,0,   &
            var(i1,0,0),1,xhalo,right,0, &
            comm_cart,status,error)
       call MPI_SENDRECV(var(imax,0,0),1,xhalo,right,0, &
            var(0,0,0),1,xhalo,left,0,     &
            comm_cart,status,error)
       !call MPI_IRECV(var(0,0,0),1,xhalo,left,1, &
       !               comm_cart,requests(2),error)
       !call MPI_IRECV(var(i1,0,0),1,xhalo,right,0, &
       !               comm_cart,requests(1),error)
       !call MPI_ISSEND(var(imax,0,0),1,xhalo,right,1, &
       !               comm_cart,requests(4),error)
       !call MPI_ISSEND(var(1,0,0),1,xhalo,left,0, &
       !               comm_cart,requests(3),error)
       !call MPI_WAITALL(4, requests, statuses, error)
    case(2) ! y direction
       call MPI_SENDRECV(var(0,1,0),1,yhalo,front,0, &
            var(0,j1,0),1,yhalo,back,0, &
            comm_cart,status,error)
       call MPI_SENDRECV(var(0,jmax,0),1,yhalo,back,0, &
            var(0,0,0),1,yhalo,front,0,   &
            comm_cart,status,error)
       !call MPI_IRECV(var(0,j1,0),1,yhalo,back,0, &
       !               comm_cart,requests(1),error)
       !call MPI_IRECV(var(0,0,0),1,yhalo,front,1, &
       !               comm_cart,requests(2),error)
       !call MPI_ISSEND(var(0,1,0),1,yhalo,front,0, &
       !               comm_cart,requests(3),error)
       !call MPI_ISSEND(var(0,jmax,0),1,yhalo,back,1, &
       !               comm_cart,requests(4),error)
       !call MPI_WAITALL(4, requests, statuses, error)
    end select
    !
    return
  end subroutine updthalos
  !
  subroutine updthalos_ibm(var,dir,icase)
    use mpi
    use mod_param
    use mod_common_mpi
    implicit none
    real, dimension(-1:,-1:,-1:), intent(inout) :: var !change later !!!
    integer, intent(in) :: dir,icase
    !real, dimension(-1:i1+1,-1:j1+1,-1:k1+1) :: hal ! THIS WILL BE CHANGED TO PLANES
    real, dimension(1:2,-1:j1+1,-1:k1+1) :: xhm
    real, dimension(imax-1:imax,-1:j1+1,-1:k1+1) :: xhp
    real, dimension(-1:i1+1,1:2,-1:k1+1) :: yhm
    real, dimension(-1:i1+1,jmax-1:jmax,-1:k1+1) :: yhp
    !integer :: requests(4), statuses(MPI_STATUS_SIZE,4)
    !
    !  This subroutine updates the halos that store info
    !  from the neighboring computational sub-domain
    !
    select case(icase)
    case(1)
       select case(dir)
       case(1) ! x direction
          call MPI_SENDRECV(var(1 ,-1,-1),1,xhalo_ibm,left ,0, &
               var(i1,-1,-1),1,xhalo_ibm,right,0, &
               comm_cart,status,error)
          call MPI_SENDRECV(var(imax-1,-1,-1),1,xhalo_ibm,right,0, &
               var(-1    ,-1,-1),1,xhalo_ibm,left ,0, &
               comm_cart,status,error)
          !    call MPI_IRECV(var(i1,-1,-1),1,xhalo_ibm,right,1, &
          !                   comm_cart,requests(1),error)
          !    call MPI_IRECV(var(-1,-1,-1),1,xhalo_ibm,left,0, &
          !                   comm_cart,requests(2),error)
          !    call MPI_ISSEND(var(imax-1,-1,-1),1,xhalo_ibm,right,0, &
          !                   comm_cart,requests(3),error)
          !    call MPI_ISSEND(var(1,-1,-1),1,xhalo_ibm,left,1, &
          !                   comm_cart,requests(4),error)
          !    call MPI_WAITALL(4, requests, statuses, error)
       case(2) ! y direction
          call MPI_SENDRECV(var(-1,1 ,-1),1,yhalo_ibm,front,0, &
               var(-1,j1,-1),1,yhalo_ibm,back ,0, &
               comm_cart,status,error)
          call MPI_SENDRECV(var(-1,jmax-1,-1),1,yhalo_ibm,back ,0, &
               var(-1,-1    ,-1),1,yhalo_ibm,front,0, &
               comm_cart,status,error)
          !    call MPI_IRECV(var(-1,j1,-1),1,yhalo_ibm,back ,1, &
          !                   comm_cart,requests(1),error)
          !    call MPI_IRECV(var(-1,-1,-1),1,yhalo_ibm,front,0, &
          !                   comm_cart,requests(2),error)
          !    call MPI_ISSEND(var(-1,jmax-1,-1),1,yhalo_ibm,back,0, &
          !                   comm_cart,requests(3),error)
          !    call MPI_ISSEND(var(-1,1,-1),1,yhalo_ibm,front,1, &
          !                   comm_cart,requests(4),error)
          !    call MPI_WAITALL(4, requests, statuses, error)
       end select
    case(2)
       select case(dir)
       case(1) ! x direction
          call MPI_SENDRECV(var(i1,-1,-1),1,xhalo_ibm,right,0, &
               xhm(1 ,-1,-1),(j1+1+2)*(k1+1+2)*2,MPI_REAL8,left ,0, &
               comm_cart,status,error)
          !$omp workshare
          var(1:2,-1:j1+1,-1:k1+1) = var(1:2,-1:j1+1,-1:k1+1) + xhm(1:2,-1:j1+1,-1:k1+1)
          !$omp end workshare
          call MPI_SENDRECV(var(-1    ,-1,-1),1,xhalo_ibm,left ,0, &
               xhp(imax-1,-1,-1),(j1+1+2)*(k1+1+2)*2,MPI_REAL8,right,0, &
               comm_cart,status,error)
          !$omp workshare
          var(imax-1:imax,-1:j1+1,-1:k1+1) = var(imax-1:imax,-1:j1+1,-1:k1+1) + xhp(imax-1:imax,-1:j1+1,-1:k1+1)
          !$omp end workshare
          !    
       case(2) ! y direction
          call MPI_SENDRECV(var(-1,j1,-1),1,yhalo_ibm,back ,0, &
               yhm(-1,1 ,-1),(i1+1+2)*(k1+1+2)*2,MPI_REAL8,front,0, &
               comm_cart,status,error)
          !$omp workshare
          var(-1:i1+1,1:2,-1:k1+1) = var(-1:i1+1,1:2,-1:k1+1) + yhm(-1:i1+1,1:2,-1:k1+1)
          !$omp end workshare
          call MPI_SENDRECV(var(-1,-1    ,-1),1,yhalo_ibm,front,0, &
               yhp(-1,jmax-1,-1),(i1+1+2)*(k1+1+2)*2,MPI_REAL8,back ,0, &
               comm_cart,status,error)
          !$omp workshare
          var(-1:i1+1,jmax-1:jmax,-1:k1+1) = var(-1:i1+1,jmax-1:jmax,-1:k1+1) + yhp(-1:i1+1,jmax-1:jmax,-1:k1+1)
          !$omp end workshare
       end select
    end select
    !
    return
  end subroutine updthalos_ibm
  !
end module mod_bound
