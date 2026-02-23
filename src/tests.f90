module mod_tests
  use mod_param
  use mod_common
  use mod_bound
  use decomp_2d 
  use mpi
  use mod_initmpi
  use mod_common_mpi 
  use mod_rk, only: rk1, rk2, rk3, rk_stage
  implicit none
  private
  public tests_scalar_bc, &
         tests_scalar_scheme_L2, &
         tests_rk_stage_equivalence, &
         tests_scalar_scheme_init, &
         tests_scalar_scheme_output, &
         tests_scalar_active_output, &
         test_sphere_int_param, & ! for debugging sphere integration issues (SS)
         test_sphere_int_second_part2
contains

  ! ------------------------------------------------------------------------- !
  ! subroutines for TEST_SCALAR_BC

  subroutine tests_scalar_bc
    implicit none
    integer :: i, j, k
    character(len=60) :: fname

    if (myid == 0) then
      write(*,*) 'Running TEST_SCALAR_BC ...'
    end if

    ! write processor topology to screen
    write(*,'(A,I4,A,I3,A,I3,A,4(I4,1X))') 'Rank=',nrank, &
        ' coorankx=',coords(1),' cooranky=',coords(2),' neighbors L/R/F/B=', left, right, front, back


    ! test x-directions
    write(fname,'("scalar_x_",I4.4,"_", I4.4,".dat")') coords(1), coords(2)

    ! Initialize the fields
    do k=0,k1
        do j = 0, j1
            do i = 0, i1
                cnew(i,j,k) =  xc(i)
            end do
        end do
    end do
    call boundscal(cnew) 

    open(unit=11, file=fname, status='replace', action='write')   
    do k = 0, k1
      do j = 0, j1
        write(11, '(1000E12.4)') (cnew(i,j,k), i=0, i1)
      end do
    end do
    close(11)

    ! test y-directions
    write(fname,'("scalar_y_",I4.4,"_", I4.4,".dat")') coords(1), coords(2)

    ! Initialize the fields
    do k=0,k1
        do j = 0, j1
            do i = 0, i1
                cnew(i,j,k) =  yc(j)
            end do
        end do
    end do
    call boundscal(cnew) 

    open(unit=11, file=fname, status='replace', action='write')   
    do k = 0, k1
      do j = 0, j1
        write(11, '(1000E12.4)') (cnew(i,j,k), i=0, i1)
      end do
    end do
    close(11)

    ! test z-directions
    write(fname,'("scalar_z_",I4.4,"_", I4.4,".dat")') coords(1), coords(2)

    ! Initialize the fields
    do k=0,k1
        do j = 0, j1
            do i = 0, i1
                cnew(i,j,k) =  zc(k)
            end do
        end do
    end do
    call boundscal(cnew) 

    open(unit=11, file=fname, status='replace', action='write')   
    do k = 0, k1
      do j = 0, j1
        write(11, '(1000E12.4)') (cnew(i,j,k), i=0, i1)
      end do
    end do
    close(11)

  end subroutine tests_scalar_bc

  ! ------------------------------------------------------------------------- !
  ! subroutines for TEST_SCALAR_SCHEME

  subroutine tests_scalar_scheme_init(c) 
    implicit none
    real, dimension(0:,0:,0:), intent(out) :: c
    integer :: i, j, k

    do k = 1,kmax
      do j = 1,jmax
        do i = 1,imax
           c(i,j,k) = sin(pi*xc(i)) * sin(pi*yc(j)) * sin(pi*zc(k))
        end do
      end do
    end do
    call boundscal(c)

  end subroutine tests_scalar_scheme_init

  ! ------------------------------------------------------------------------- !

  ! tests_scalar_scheme_src moved to mod_test_sources to avoid module cycle

  subroutine tests_scalar_scheme_output(unew, vnew, wnew, cnew,pnew, istep)
    implicit none
    real, dimension(0:,0:,0:), intent(in) :: unew, vnew, wnew, cnew,pnew
    integer :: i, j, k, istep
    character(len=60) :: fname
    write(fname,'("unew_rank","_",I4.4,"_",I4.4,"_",I4.4,".dat")') coords(1), coords(2), istep
    open(unit=1000, file=fname, status='replace', action='write')
    do k = 0, k1
      do j = 0, j1
        write(1000, '(1000E12.4)') (unew(i,j,k), i=0, i1)
      end do
    end do 
    close(1000)

    write(fname,'("vnew_rank","_",I4.4,"_",I4.4,"_",I4.4,".dat")') coords(1), coords(2), istep
    open(unit=1001, file=fname, status='replace', action='write')
    do k = 0, k1
      do j = 0, j1
        write(1001, '(1000E12.4)') (vnew(i,j,k), i=0, i1)
      end do
    end do 
    close(1001)

    write(fname,'("wnew_rank","_",I4.4,"_",I4.4,"_",I4.4,".dat")') coords(1), coords(2), istep
    open(unit=1002, file=fname, status='replace', action='write')
    do k = 0, k1
      do j = 0, j1
        write(1002, '(1000E12.4)') (wnew(i,j,k), i=0, i1)
     end do
    end do 
    close(1002)  

    write(fname,'("cnew_rank","_",I4.4,"_",I4.4,"_",I4.4,".dat")') coords(1), coords(2), istep
    open(unit=1003, file=fname, status='replace', action='write')
    do k = 0, k1
      do j = 0, j1
        write(1003, '(1000E12.4)') (cnew(i,j,k), i=0, i1)
     end do
    end do 
    close(1003)      
  end subroutine tests_scalar_scheme_output

  ! ------------------------------------------------------------------------- !

  subroutine tests_scalar_scheme_L2(c, time)
    real, dimension(0:,0:,0:), intent(in) :: c
    real, intent(in) :: time

    real, dimension(0:i1,0:j1,0:k1) :: c_exact
    real :: local_err, global_err, L2err
    character(len=30) :: fname
    real :: coorx, coory, coorz
    integer :: i, j, k, ierror,iglob,jglob,kglob

    do k = 1,kmax
      coorz = zc(k)/lz ! normalised with channel height using grid arrays
      do j = 1,jmax
        do i = 1,imax
           coorx = xc(i)/lx
           coory = yc(j)/ly
           ! Analytical profile for current MMS test (no explicit time dependence):
           c_exact(i,j,k) = sin(pi*coorx) * sin(pi*coory) * sin(pi*coorz)
        end do
      end do
    end do

    local_err = 0.0
    do k = 1, kmax
      do j = 1, jmax
        do i = 1, imax    
          local_err = local_err + (cnew(i,j,k) - c_exact(i,j,k))**2 ! For \frac{\partial \c}{\partial t} = 0 version,
          ! cnew = c_exact always, they don't change with time. For \frac{\partial \c}{\partial t} \= 0 version, 
          ! cnew still doesn't change with time
        end do
      end do
    end do  

    call MPI_Allreduce(local_err, global_err, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
    !print *, 'local_err (rank', nrank, ') =', local_err
    !print *, 'global_err (rank', nrank, ') =', global_err
    
    L2err = sqrt(global_err / (itot*jtot*ktot))
    ! L2err = sqrt(global_err)
    if (nrank == 0) print *, 'MMS L2 Error =', L2err

    write(fname,'("analytical_rank","_",I4.4,"_",I4.4,".dat")') coords(1), coords(2)
    open(unit=20, file=fname, status='replace', action='write')
    do k = 1, kmax
      do j = 1, jmax
        write(20, '(1000E12.4)') (c_exact(i,j,k), i=1, imax)
      end do
    end do
    close(20)

  end subroutine tests_scalar_scheme_L2

  ! ------------------------------------------------------------------------- !
  ! subroutines for TEST_SCALAR_ACTIVE

  subroutine tests_scalar_active_output(unew, vnew, wnew, cnew,pnew, istep)
    implicit none
    real, dimension(0:,0:,0:), intent(in) :: unew, vnew, wnew, cnew,pnew
    integer :: i, j, k, istep, imax, jmax, kmax
    character(len=60) :: fname

  if ((nrank == 0 .or. nrank == 2) .and. mod(istep,20) == 0) then
      write(fname,'("unew_rank",I4.4,"_",I4.4,".dat")') nrank, istep
      open(unit=1200, file=fname, status='replace', action='write')
      j = (j1-1)/2 ! output middle slice of a rank in y direction
      do k = 1, k1 - 1
          write(1200, '(1000E12.4)') (unew(i,j,k), i=1, i1-1)
      end do
      close(1200)

      write(fname,'("vnew_rank",I4.4,"_",I4.4,".dat")') nrank, istep
      open(unit=1201, file=fname, status='replace', action='write')
      j = (j1-1)/2 ! output middle slice of a rank in y direction
      do k = 1, k1 - 1
          write(1201, '(1000E12.4)') (vnew(i,j,k), i=1, i1-1)
      end do
      close(1201)

      write(fname,'("wnew_rank",I4.4,"_",I4.4,".dat")') nrank, istep
      open(unit=1202, file=fname, status='replace', action='write')
      j = (j1-1)/2 ! output middle slice of a rank in y direction
      do k = 1, k1 - 1  
          write(1202, '(1000E12.4)') (wnew(i,j,k), i=1, i1-1)
      end do
      close(1202)

      write(fname,'("cnew_rank",I4.4,"_",I4.4,".dat")') nrank, istep
      open(unit=1203, file=fname, status='replace', action='write')
      j = (j1-1)/2 ! output middle slice of a rank in y direction
      do k = 1, k1 - 1
          write(1203, '(1000E12.4)') (cnew(i,j,k), i=1, i1-1)
      end do
      close(1203)

      if (myid .eq. 0) then
      open(unit=1204, file='RBC_time_record.txt', status='unknown', action='write', position='append')
      write(1204,*) time
      close(1204)
      end if
  end if

  if (mod(istep,20) == 0) then
      write(fname,'("cnew_thermal_layer_rank",I4.4,"_",I4.4,".dat")') nrank, istep
      open(unit=1205, file=fname, status='replace', action='write')
      k = 1 ! output thermal layer close to bottom wall
      do j = 1, j1 - 1
          write(1205, '(1000E12.4)') (cnew(i,j,k), i=1, i1-1)
      end do
      close(1205)

      write(fname,'("cnew_center_rank",I4.4,"_",I4.4,".dat")') nrank, istep
      open(unit=1206, file=fname, status='replace', action='write')
      k = (k1-1)/2 ! output middle slice in z direction
      do j = 1, j1 - 1
          write(1206, '(1000E12.4)') (cnew(i,j,k), i=1, i1-1)
      end do
      close(1206)

  end if


  end subroutine tests_scalar_active_output

  ! ------------------------------------------------------------------------- !
  ! Test: verify three rk_stage passes are bitwise-identical to rk1->rk2->rk3
  ! To avoid a module cycle (mod_rk uses mod_tests), this test accepts
  ! two runner procedures provided by the caller:
  !  - run_classic: must call rk1, rk2, rk3 on the given arrays
  !  - run_stage:   must call rk_stage(0), rk_stage(1), rk_stage(2)
  ! The test sets identical initial conditions and compares dudt,dvdt,dwdt,dcdt
  ! and wall-shear bookkeeping exactly (bitwise via .eq.).

  ! Shuang: it is bitwise identical. 
  subroutine tests_rk_stage_equivalence()
    implicit none
    ! locals
    real, dimension(0:i1,0:j1,0:k1) :: durk, dvrk, dwrk, dcrk
    real, dimension(0:i1,0:j1,0:k1) :: dudt_ref, dvdt_ref, dwdt_ref, dcdt_ref

    real :: ws_old_init, ws_new_init, ws_old_ref, ws_new_ref
    integer :: i,j,k, ierror,ii
    integer(kind=8) :: mism_local, mism_global
    character(len=64) :: fname1101, fname1102

    ! deterministic initial conditions (velocity: TG vortex at t=0; scalar: sin-sin-sin)
    call tests_scalar_scheme_init(cnew)
      pnew(:,:,:) = 0.0
    ! time = 0.0
    ! dt   = 1.0e-3

    durk = 0.0; dvrk = 0.0; dwrk = 0.0; dcrk = 0.0
    dudt = 0.0; dvdt = 0.0; dwdt = 0.0; dcdt = 0.0
    ws_old_init = wallshearold
    ws_new_init = wallshearnew
    wallshearold = 0.0
    wallshearnew = 0.0

  ! path A: classic rk1->rk2->rk3
  call rk1(durk,dvrk,dwrk,dcrk) 
  call rk2(durk,dvrk,dwrk,dcrk)
  call rk3(durk,dvrk,dwrk,dcrk)

    dudt_ref = dudt; dvdt_ref = dvdt; dwdt_ref = dwdt; dcdt_ref = dcdt
    ws_old_ref = wallshearold
    ws_new_ref = wallshearnew

    ! reset state
    durk = 0.0; dvrk = 0.0; dwrk = 0.0; dcrk = 0.0
    dudt = 0.0; dvdt = 0.0; dwdt = 0.0; dcdt = 0.0
    wallshearold = 0.0
    wallshearnew = 0.0

  ! path B: rk_stage(0/1/2)
  call rk_stage(0,durk,dvrk,dwrk,dcrk)
  call rk_stage(1,durk,dvrk,dwrk,dcrk)
  call rk_stage(2,durk,dvrk,dwrk,dcrk)


    ! compare bitwise
    mism_local = 0
    do k=0,k1
      do j=0,j1
        do i=0,i1
          if (dudt(i,j,k) .ne. dudt_ref(i,j,k)) mism_local = mism_local + 1
          if (dvdt(i,j,k) .ne. dvdt_ref(i,j,k)) mism_local = mism_local + 1
          if (dwdt(i,j,k) .ne. dwdt_ref(i,j,k)) mism_local = mism_local + 1
          if (dcdt(i,j,k) .ne. dcdt_ref(i,j,k)) mism_local = mism_local + 1
        end do
      end do
    end do
    if (wallshearold .ne. ws_old_ref) mism_local = mism_local + 1
    if (wallshearnew .ne. ws_new_ref) mism_local = mism_local + 1


    call MPI_Allreduce(mism_local, mism_global, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierror)

    if (nrank == 0) then
      if (mism_global == 0) then
        write(*,*) 'RK stage equivalence: PASS (bitwise identical).'
      else
        write(*,*) 'RK stage equivalence: FAIL; mismatches = ', mism_global
      end if
    end if

  ! Test dcdt_ref and dcdt
  ! write(fname1101,'("classic_dcdt",I4.4,"_",I4.4,".dat")') coords(1), coords(2)
  ! open(unit=1101, file=fname1101, status='replace', action='write')
  ! do k = 0, k1
  !    do j = 0, j1
  !       write(1101, '(1000E12.4)') (dcdt_ref(i,j,k), i=0, i1)
  !    end do
  ! end do 
  ! close(1101)

  ! write(fname1102,'("unify_dcdt",I4.4,"_",I4.4,".dat")') coords(1), coords(2)
  ! open(unit=1102, file=fname1102, status='replace', action='write')
  ! do k = 0, k1
  !    do j = 0, j1
  !       write(1102, '(1000E12.4)') (dcdt(i,j,k), i=0, i1)
  !    end do
  ! end do 
  ! close(1102)

    ! restore inputs
    wallshearold = ws_old_init
    wallshearnew = ws_new_init
  end subroutine tests_rk_stage_equivalence

  ! ------------------------------------------------------------------------- !
  ! Debug sphere integration--check input parameters
  ! for debugging sphere integration issues (SS)
  subroutine test_sphere_int_param(radin, radin2, dx, dy, dz, type)
    implicit none
    integer, intent(in) :: type
    real, intent(in) :: radin, radin2
    real, intent(in) :: dx,dy,dz
    ! Debug output (SS_debug)
    if (myid == 10 .and. type == 1) then
       write(6,*) ''
       write(6,*) '=== INT SPHERE PARAMETER TEST ==='
       write(6,*) 'radius = ', radius
       write(6,*) 'dx,dy,dz = ', dx, dy, dz
       write(6,*) 'radin = ', radin, ' radin2 = ', radin2
       write(6,*) 'dVeul = ', dVeul
       write(6,*) 'pmax = ', pmax
       write(6,*) '=== END TEST ==='
       write(6,*) ''
    endif
    
  end subroutine test_sphere_int_param

  ! ------------------------------------------------------------------------- !
  ! Debug sphere integration--check first step and second step part 1 results
  ! Put these two tests in the intgr_over_sphere.f90 because the anb have three different definitions
  ! in forcing.f90, phase_indicator.f90 and intgr_over_sphere.f90.
  ! for debugging sphere integration issues (SS)  

  ! ------------------------------------------------------------------------- !
  ! Debug sphere integration--check second step part 2 results
  ! for debugging sphere integration issues (SS)
  subroutine test_sphere_int_second_part2(nb, p, coorxc,cooryc,coorzc,ilow,ihigh,jlow,jhigh,klow,khigh,type)
    implicit none
    integer, intent(in) :: type
    integer, intent(in) :: nb, p
    integer, intent(in) :: ilow,ihigh,jlow,jhigh,klow,khigh
    real, intent(in) :: coorxc,cooryc,coorzc    
    if (myid == 10 .and. type == 1 .and. p == 1 .and. nb == 0) then
       write(6,*) ''
       write(6,*) '=== INT SPHERE SECOND STEP PART 2 TEST ==='
       open(unit=99, file='debug_integration.txt', status='replace')
       write(99,'(A)') 'nb     coorxc     cooryc     coorzc   ilow ihigh jlow jhigh klow khigh'
       write(99,'(A)') '---------------------------------------------------------------------------------'
       write(99,'(I2, 3F12.4, 6I6)') nb,  coorxc, cooryc, coorzc, ilow, ihigh, jlow, jhigh, klow, khigh
       close(99)
    
       write(6,*) '=== END SECOND STEP PART 2 TEST ==='
       write(6,*) ''
    end if

  end subroutine test_sphere_int_second_part2  

end module mod_tests
