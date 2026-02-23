!                                                                                                                                                                                        
!                                                                                                                                                                                        
! IIIIIIIIII                         tttt                            RRRRRRRRRRRRRRRRR   PPPPPPPPPPPPPPPPP                                               tttt             SSSSSSSSSSSSSSS 
! I::::::::I                      ttt:::t                            R::::::::::::::::R  P::::::::::::::::P                                           ttt:::t           SS:::::::::::::::S
! I::::::::I                      t:::::t                            R::::::RRRRRR:::::R P::::::PPPPPP:::::P                                          t:::::t          S:::::SSSSSS::::::S
! II::::::II                      t:::::t                            RR:::::R     R:::::RPP:::::P     P:::::P                                         t:::::t          S:::::S     SSSSSSS
!   I::::Innnn  nnnnnnnn    ttttttt:::::ttttttt        eeeeeeeeeeee    R::::R     R:::::R  P::::P     P:::::Paaaaaaaaaaaaa  rrrrr   rrrrrrrrr   ttttttt:::::ttttttt    S:::::S            
!   I::::In:::nn::::::::nn  t:::::::::::::::::t      ee::::::::::::ee  R::::R     R:::::R  P::::P     P:::::Pa::::::::::::a r::::rrr:::::::::r  t:::::::::::::::::t    S:::::S            
!   I::::In::::::::::::::nn t:::::::::::::::::t     e::::::eeeee:::::eeR::::RRRRRR:::::R   P::::PPPPPP:::::P aaaaaaaaa:::::ar:::::::::::::::::r t:::::::::::::::::t     S::::SSSS         
!   I::::Inn:::::::::::::::ntttttt:::::::tttttt    e::::::e     e:::::eR:::::::::::::RR    P:::::::::::::PP           a::::arr::::::rrrrr::::::rtttttt:::::::tttttt      SS::::::SSSSS    
!   I::::I  n:::::nnnn:::::n      t:::::t          e:::::::eeeee::::::eR::::RRRRRR:::::R   P::::PPPPPPPPP      aaaaaaa:::::a r:::::r     r:::::r      t:::::t              SSS::::::::SS  
!   I::::I  n::::n    n::::n      t:::::t          e:::::::::::::::::e R::::R     R:::::R  P::::P            aa::::::::::::a r:::::r     rrrrrrr      t:::::t                 SSSSSS::::S 
!   I::::I  n::::n    n::::n      t:::::t          e::::::eeeeeeeeeee  R::::R     R:::::R  P::::P           a::::aaaa::::::a r:::::r                  t:::::t                      S:::::S
!   I::::I  n::::n    n::::n      t:::::t    tttttte:::::::e           R::::R     R:::::R  P::::P          a::::a    a:::::a r:::::r                  t:::::t    tttttt            S:::::S
! II::::::IIn::::n    n::::n      t::::::tttt:::::te::::::::e        RR:::::R     R:::::RPP::::::PP        a::::a    a:::::a r:::::r                  t::::::tttt:::::tSSSSSSS     S:::::S
! I::::::::In::::n    n::::n      tt::::::::::::::t e::::::::eeeeeeeeR::::::R     R:::::RP::::::::P        a:::::aaaa::::::a r:::::r                  tt::::::::::::::tS::::::SSSSSS:::::S
! I::::::::In::::n    n::::n        tt:::::::::::tt  ee:::::::::::::eR::::::R     R:::::RP::::::::P         a::::::::::aa:::ar:::::r                    tt:::::::::::ttS:::::::::::::::SS 
! IIIIIIIIIInnnnnn    nnnnnn          ttttttttttt      eeeeeeeeeeeeeeRRRRRRRR     RRRRRRRPPPPPPPPPP          aaaaaaaaaa  aaaarrrrrrr                      ttttttttttt   SSSSSSSSSSSSSSS 
! 
! InteRPartS - Interface-Resolved Particle Simulations
! 
! Code for channel (full or half) transport
!
! Contributors: Pedro Costa & Wim-Paul Breugem
! p.simoes.costa@gmail.com / w.p.breugem@tudelft.nl
!
! Last modified: July 2015
!
!**********************************************************
program interpt
  use decomp_2d_io
  use decomp_2d
  use mod_param
  use mod_common
  use mod_common_mpi
  use mod_initmpi
  use mod_init
  use mod_bound
  use mod_chkdiv
  use mod_chkdt
  use mod_loadflds
  use mod_loadpart
  use mod_rk
  use mod_initsolver
  use mod_fillps
  use mod_solver
  use mod_correc
  use mod_coordsfp
  use mod_initparticles
  use mod_kernel
  use mod_intgr_over_sphere
  use mod_interp_spread
  use mod_forcing
  use mod_scal ! (SS)
  use mod_tests
  use mod_intgr_nwtn_eulr
  use mod_intgr_nwtn_eulr_subc
  use mod_output
  use mod_fftw
  use mod_autocorr
  use mod_tests

  implicit none
  integer :: iopt,begin,nstep,istep,nsub
  real :: dtmax
  real, dimension(0:i1,0:j1,0:k1) :: durkold,dvrkold,dwrkold,dcrkold
  real, target, dimension(0:i1,0:j1,0:k1) :: p
  real, pointer, dimension(:,:,:) :: pp ! uncomment iff solver2d is used
  real :: maxerror,avererror,maxerrorsum,avererrorsum
  integer :: ibmiter
  integer :: i,j,k
  integer :: ierror
  real :: norm,norm_all,sumk
  real :: xpos, ypos, zpos
  !real, parameter :: pi = 3.141592653589793
  real :: v_bulk_all
  character(len=4) :: number
  real :: comp_time,comp_time_avrg,comp_time_min,comp_time_max ! for mpi timer
  real :: comp_time_fluid,comp_time_ibm,comp_time_nwtn_eulr,comp_time_aux
  real :: tsave,ttemp,tmax
  integer :: ihistory,ihistory_acc,count_history,count_history_max
  integer :: iamplify
  real :: factor
  real :: dpdy

  !real, allocatable, dimension(:,:,:) :: temp (COMMENTED OUT BY Pedro: COBP)
  !
  iopt = 0
  if(iopt.eq.1) then
     call MPI_INIT(error)
     call decomp_2d_init(itot, jtot, ktot, 0,0)
     call MPI_FINALIZE(error)
     stop
  else
     call initmpi
  endif

  call input_parameters
  call init_grid
  ! call derived_parameters
! Output radius parameter
if (np > 0) then  
if (myid .eq. 0) then
   ! write(6,*) 'Particle radius = ', radius
   ! write(6,*) 'Kinematic viscosity = ', visc
   ! write(6,*) 'post1d frequency = ', iout1d
   ! write(6,*) 'Gravity acceleration components:'
   ! write(6,*) '  gaccx = ', gaccx
   ! write(6,*) '  gaccy = ', gaccy
   ! write(6,*) '  gaccz = ', gaccz
   write(6,*) 'Number of Lagrangian Force Points of first shell = ', nl
   write(6,*) 'Number of Lagrangian Force Points of second shell = ', NL2
   write(6,*) 'Number of Lagrangian Force Points of third shell = ', NL3
   write(6,*) 'Number of Lagrangian Force Points of fourth shell = ', NL4
   write(6,*) 'Number of Lagrangian Force Points of all shells = ', NLtot   
endif
! Check if nl equals nlmax, stop simulation if not
if (nl .ne. nlmax) then
   if (myid .eq. 0) then
      write(6,*) 'ERROR: Number of Lagrangian Force Points (nl =', nl, ') does not match expected value (nlmax =', nlmax, ')'
      write(6,*) 'Simulation stopped due to mismatch in Lagrangian Force Points'
   endif
   ! optimize program stop implementation
   call decomp_2d_finalize ! Clean up: release memory and resources allocated by 2DECOMP&FFT library
   call MPI_FINALIZE(error) ! finalize MPI
   ! call MPI_ABORT(comm_cart, error_code, error) ! force all processes to terminate.
   stop
else
   if (myid .eq. 0) then
      write(6,*) 'SUCCESS: nl matches nlmax =', nlmax
   endif
endif
endif
  ! write(*,*) myid, runmode
  !
  tsave = MPI_WTIME()       
  begin = 0        
  nstep = 100000             
  nsub = 200
  tmax = 23.5 ! hours
  if(begin.eq.0) then  ! THIS IS ACTUALLY A SWITCH FOR COLD/WARM START
     call init_ics         
     if (np > 0) then
        if (myid .eq. 0) then
            call kerneltest(sumk)
            write(6,*) 'Integral over kernel = ', sumk
        endif  
        call initparticles
     end if 
     time = 0.
     call bounduvw(unew,vnew,wnew)
     call boundscal(cnew)
     call chkdt(dtmax)
     dt = 0.5*dtmax
  else
     !$omp workshare
     unew(:,:,:) = 0.
     vnew(:,:,:) = 0.
     wnew(:,:,:) = 0.
     pnew(:,:,:) = 0.
     cnew(:,:,:) = 0.
     !$omp end workshare
     if (np > 0) then  ! Restore the calling order of Pedro's original files (Shuang Shan: SS)
        call loadpart(0,begin)
     end if
     call loadflds(0,begin) ! 0:read, 1:write (NOTE: this will overwrite begin with value from file) (COBP)
     !  allocate(temp(1:imax,1:jmax,1:kmax)) ! (COBP) in the next 10 lines
     !  temp(:,:,:) = unew(1:imax,1:jmax,1:kmax)
     !  call decomp_2d_write_every(3,temp,2,2,2,'unew_smaller.bin',.true.)
     !  temp(:,:,:) = vnew(1:imax,1:jmax,1:kmax)
     !  call decomp_2d_write_every(3,temp,2,2,2,'vnew_smaller.bin',.true.)
     !  temp(:,:,:) = wnew(1:imax,1:jmax,1:kmax)
     !  call decomp_2d_write_every(3,temp,2,2,2,'wnew_smaller.bin',.true.)
     !  temp(:,:,:) = pnew(1:imax,1:jmax,1:kmax)
     !  call decomp_2d_write_every(3,temp,2,2,2,'pnew_smaller.bin',.true.)
     !  deallocate(temp)
     !  call mpi_finalize(error)
     iamplify = 0
     if(iamplify.eq.1) then
        factor = 10.
        call amplify_vel_fluctuations(factor)
        call bounduvw(unew,vnew,wnew)
        call boundscal(cnew)
        call chkdt(dtmax)
        dt = min(dt,0.5*dtmax)
     endif
      call chkdt(dtmax) ! (COBP) in the next 3 lines
      ! dt = 0.5*dtmax
      ! dt = 0.75*dtmax
      ! dt = min(0.050,dt)
     if (myid .eq. 0) write(6,*) 'nr steps at beginning simulation = ',begin
  endif
  !
  if (myid .eq. 0) write(6,*) 'Solid volume fraction of particles = ',solidity
  call initsolver
  pp => p(1:imax,1:jmax,1:kmax) ! uncomment iff solver2d is used
  call bounduvw(unew,vnew,wnew)
  call boundscal(cnew)
  call boundp(pnew)
  call chkdiv 
  if (np > 0) then
     call coordsfp  
     if(begin.eq.0) then
        call intgr_over_sphere(1) 
        call intgr_over_sphere(2)
        call intgr_over_sphere(3)
        !call intgr_over_sphere(ap(:)%intu,ap(:)%intv,ap(:)%intw,1) ! (COBP) in the next 2 lines
        !call intgr_over_sphere(ap(:)%intu,ap(:)%intv,ap(:)%intw,2)
        !call intgr_over_sphere(ap(:)%intomx,ap(:)%intomy,ap(:)%intomz,3)
     endif
  end if

  !
  ! main time-integration loop below
  !
  
  if (myid .eq. 0) write(6,*) 'dtmax = ', dtmax, ' dt = ', dt
  !$omp workshare
  dudt(:,:,:) = 0.
  dvdt(:,:,:) = 0.
  dwdt(:,:,:) = 0.
  dcdt(:,:,:) = 0.
  durkold(:,:,:) = 0.
  dvrkold(:,:,:) = 0.
  dwrkold(:,:,:) = 0.
  dcrkold(:,:,:) = 0.
  !$omp end workshare
  !
  ! initialize fftw for computing autocorrelations
  ! Restore the calling order of Pedro's original files (SS)
  call init_fft(itot,plan_r2c_x,plan_c2r_x)
  call init_fft(jtot,plan_r2c_y,plan_c2r_y)
  ihistory = 0
  ihistory_acc = 0

  ! switch for simple unit tests
  select case (runmode)
  case (RUNMODE_TEST_RK_STAGE)
     call tests_rk_stage_equivalence()
     call decomp_2d_finalize
     call MPI_FINALIZE(error)
     stop   
  case (RUNMODE_TEST_SCALAR_BCS)
     call tests_scalar_bc()
     call decomp_2d_finalize
     call MPI_FINALIZE(error)
     stop
  case (RUNMODE_TEST_SCALAR_SCHEME)
     ! write output to check initial condition
     call tests_scalar_scheme_output(unew, vnew, wnew, cnew,pnew, 0)
   end select
  
  do istep = begin+1,nstep
     if (time + dt .gt. t_end) then
        dt = t_end - time
     endif
     time = time + dt
     if (myid.eq.0) write(6,*) 'time = ', time, 'istep = ', istep
     maxerrorsum  = 0.
     avererrorsum = 0.
     comp_time = MPI_WTIME()
     !$omp workshare
     v_bulk =sum(vnew(1:imax,1:jmax,1:kmax))
     !$omp end workshare
     call mpi_allreduce(v_bulk,v_bulk_all,1,mpi_real8,mpi_sum,comm_cart,error)
     v_bulk=v_bulk_all/(1.*itot*jtot*kmax)
     !     if(myid.eq.0) print*, 'Bulk velocity: ', v_bulk ! (COBP)
     do rkiter = 1,3
        rkparalpha = rkcoeffab(rkiter)
        !    comp_time = MPI_WTIME() ! (COBP)
        comp_time_fluid = MPI_WTIME()
        if (rkiter .eq. 1) then
           dpdy = 0.
           !if(myid.eq.0) write(6,*) 'RK3 step 1' ! (SS)
           call rk1(durkold,dvrkold,dwrkold,dcrkold)
        endif
        if (rkiter .eq. 2) then
           !if(myid.eq.0) write(6,*) 'RK3 step 2' ! (SS)  
           call rk2(durkold,dvrkold,dwrkold,dcrkold)
        endif
        if (rkiter .eq. 3) then
           !if(myid.eq.0) write(6,*) 'RK3 step 3' ! (SS)
           call rk3(durkold,dvrkold,dwrkold,dcrkold)
        endif

      if (np > 0) then
        !$omp workshare
        dudtold(:,:,:) = dudt(:,:,:) ! forcing.f90, interp_spread.f90: lag forces
        dvdtold(:,:,:) = dvdt(:,:,:)
        dwdtold(:,:,:) = dwdt(:,:,:)
        ! dcdtold(:,:,:) = dcdt(:,:,:)
        !$omp end workshare
      end if  
      comp_time_fluid = MPI_WTIME() - comp_time_fluid ! time used for calculating fluid parts (SS)

      if (np > 0) then
        !
        comp_time_IBM = MPI_WTIME()
        ibmiter = 0
        !$omp workshare
        dudtf(1:imax,1:jmax,1:kmax) = dudt(1:imax,1:jmax,1:kmax) ! interp_spread.f90: lag forces
        dvdtf(1:imax,1:jmax,1:kmax) = dvdt(1:imax,1:jmax,1:kmax)
        dwdtf(1:imax,1:jmax,1:kmax) = dwdt(1:imax,1:jmax,1:kmax)
        ! dcdtf(1:imax,1:jmax,1:kmax) = dcdt(1:imax,1:jmax,1:kmax) ! Scalar has not yet been 
        ! added to the fluid-structure interaction part, so it is commented out for now (SS).
        !$omp end workshare
        do while(ibmiter.le.2)
           call eulr2lagr(ibmiter)
           if (ibmiter .eq. 0) then
              call complagrforces
           else
              call updtlagrforces(maxerror,avererror)
           endif
           call lagr2eulr
           call updtintermediatevel(istep) ! (COBP)
           ibmiter = ibmiter+1
           !$omp workshare
       !    dudtf(1:imax,1:jmax,1:kmax) = dudtf(1:imax,1:jmax,1:kmax) + dudtold(1:imax,1:jmax,1:kmax) ! (COBP) in the next 2 lines
       !    dvdtf(1:imax,1:jmax,1:kmax) = dvdtf(1:imax,1:jmax,1:kmax) + dvdtold(1:imax,1:jmax,1:kmax)
       !    dwdtf(1:imax,1:jmax,1:kmax) = dwdtf(1:imax,1:jmax,1:kmax) + dwdtold(1:imax,1:jmax,1:kmax)
           !$omp end workshare
           !!!$omp workshare ! (COBP) in the next 7 lines
           !!v_bulk =sum(dvdtf(1:imax,1:jmax,1:kmax))
           !!!$omp end workshare
           !!call mpi_allreduce(v_bulk,v_bulk_all,1,mpi_real8,mpi_sum,comm_cart,error)
           !!v_bulk=v_bulk_all/(1.*itot*jtot*kmax)
           !!!$omp workshare
           !!dvdtf(1:imax,1:jmax,1:kmax) = dvdtf(1:imax,1:jmax,1:kmax) + 1.0*(bulk_v_sup-v_bulk)
           !!!$omp end workshare
        enddo
        !$omp workshare
        dudt(1:imax,1:jmax,1:kmax) = dudtf(1:imax,1:jmax,1:kmax)
        dvdt(1:imax,1:jmax,1:kmax) = dvdtf(1:imax,1:jmax,1:kmax)
        dwdt(1:imax,1:jmax,1:kmax) = dwdtf(1:imax,1:jmax,1:kmax)
        ! dcdt(1:imax,1:jmax,1:kmax) = dcdtf(1:imax,1:jmax,1:kmax)
        ! Scalar has not yet been added to the fluid-structure interaction part, 
        ! so it is commented out for now (SS).
        !$omp end workshare
        comp_time_IBM = MPI_WTIME() - comp_time_IBM
        !    ! (COBP) in the next 12 lines
        !    if (rkiter .eq. 3.and.myid.eq.0.and.np.ne.0) then
        !    maxerrorsum = maxerrorsum + rkcoeffab(rkiter)*maxerror
        !    avererrorsum = avererrorsum + rkcoeffab(rkiter)*avererror
        !      write(6,'(A25,I2,3E16.8)') 'ibmiter,max,aver,ratio = ', &
        !                                 ibmiter-1,maxerrorsum,avererrorsum,maxerrorsum/avererrorsum
        !      if (mod(istep,ioutchk).eq.0) then
        !        open(22,file=datadir//'ibmerror.txt',position='append')
        !          write(22,'(I2,4E16.8)') ibmiter-1,time/tref,maxerrorsum,avererrorsum,maxerrorsum/avererrorsum
        !        close(22)
        !      endif
        !    endif
        !
        comp_time_aux = MPI_WTIME()
      end if

      call bounduvw(dudt,dvdt,dwdt)
      call fillps(rkiter,p)
      call solver2d(pp) ! pp => p in line 160
      !    call solver1d(p)  ! (COBP)
      call boundp(p)
      call correc(rkiter,p) ! u,v,w pressure-corrected inside routine

      ! ! if using the 1001 case, update the velocity field to the exact solution at current time
      ! if (runmode == RUNMODE_TEST_SCALAR_SCHEME) then
      !   call tests_scalar_scheme_uvw(unew,vnew,wnew, time) ! only question is which time.
      ! end if

      call bounduvw(unew,vnew,wnew)
      !$omp workshare
      pnew(:,:,:) = pnew(:,:,:) + p(:,:,:)
      !$omp end workshare
      call boundp(pnew)

      ! update scalar field
      cnew(:,:,:) = dcdt(:,:,:)
      call boundscal(cnew)
    
!      k=1 ! near bottom wall  ! (Maarten van Reeuwijk: MvR) in the next 8 lines
!      !$omp workshare
!      norm = sum(sum(pnew(:,:,k),2),1)
!      !$omp end workshare
!      call mpi_allreduce(norm,norm_all,1,mpi_real8,mpi_sum,comm_cart,error)
!      norm = norm_all/(1.*itot*jtot)
!      !$omp workshare
!      pnew(:,:,:) = pnew(:,:,:) - norm
!      !$omp end workshare

      if (np > 0) then
        comp_time_fluid = comp_time_fluid + MPI_WTIME()-comp_time_aux
        comp_time_nwtn_eulr = MPI_WTIME()
        call intgr_nwtn_eulr
        !    do k=1,nsub ! (COBP) in the next 2 lines
        !      call intgr_nwtn_eulr_subc(k,nsub,dt/nsub)
        !    enddo
        comp_time_nwtn_eulr = MPI_WTIME() - comp_time_nwtn_eulr
        call mpi_allreduce(comp_time_fluid,comp_time_avrg,1,mpi_real8,mpi_sum,comm_cart,error)
        call mpi_allreduce(comp_time_fluid,comp_time_min ,1,mpi_real8,mpi_min,comm_cart,error)
        call mpi_allreduce(comp_time_fluid,comp_time_max ,1,mpi_real8,mpi_max,comm_cart,error)
        if(myid.eq.0) write(6,'(A,3E16.8)') 'Avrg, min & max elapsed time (fluid) = ', &
             comp_time_avrg/product(dims),comp_time_min,comp_time_max
        call mpi_allreduce(comp_time_IBM,comp_time_avrg,1,mpi_real8,mpi_sum,comm_cart,error)
        call mpi_allreduce(comp_time_IBM,comp_time_min ,1,mpi_real8,mpi_min,comm_cart,error)
        call mpi_allreduce(comp_time_IBM,comp_time_max ,1,mpi_real8,mpi_max,comm_cart,error)
        if(myid.eq.0) write(6,'(A,3E16.8)') 'Avrg, min & max elapsed time = (ibm)', &
             comp_time_avrg/product(dims),comp_time_min,comp_time_max
        call mpi_allreduce(comp_time_nwtn_eulr,comp_time_avrg,1,mpi_real8,mpi_sum,comm_cart,error)
        call mpi_allreduce(comp_time_nwtn_eulr,comp_time_min ,1,mpi_real8,mpi_min,comm_cart,error)
        call mpi_allreduce(comp_time_nwtn_eulr,comp_time_max ,1,mpi_real8,mpi_max,comm_cart,error)
        if(myid.eq.0) write(6,'(A,3E16.8)') 'Avrg, min & max elapsed time = (n-e)', &
             comp_time_avrg/product(dims),comp_time_min,comp_time_max
        !forceytot = -1.0*(bulk_v_sup-v_bulk)/dt ! (COBP)
        dpdy = dpdy ! + forceytot
      end if
     enddo

   !   if (runmode == RUNMODE_TEST_SCALAR_SCHEME) then
   !     call tests_scalar_scheme_l2(cnew, time) ! Calculate L2 error (test/1001: MMS)
   !     call tests_scalar_scheme_output(unew, vnew, wnew, cnew,pnew, istep)
   !   end if     

     if (runmode == RUNMODE_TEST_SCALAR_ACTIVE) then ! Rayleigh Benard Convection
      call tests_scalar_active_output(unew, vnew, wnew, cnew,pnew, istep)
     end if
       !
     if (mod(istep,ioutchk).eq.0) then
        call chkdiv
        call chkdt(dtmax)       
        !    forceytot = -1.0*(bulk_v_sup-v_bulk)/dt ! (COBP)       
        if (myid .eq. 0) then
           write(6,*) 'dtmax = ', dtmax, ' dt = ', dt
           open(29,file='timestep.txt',position='append')
           write(29,'(3E16.8)') 1.*istep,time,dt!,forceytot
           close(29)
        endif
        dt = 0.75*dtmax          
        if(myid.eq.0) then
           open(29,file='bulkforcing.txt',position='append')
           write(29,'(3E16.8)') 1.*istep,time,dpdy
           close(29)
        endif
     endif
     
     if (mod(istep,iout1d).eq.0) then
        !    do k=1,pmax  ! (COBP) in the next 10 lines
        !      if (ap(k)%mslv .gt. 0) then
        !        write(number,'(i4.4)') ap(k)%mslv
        !        open(20,file=datadir//'pchar'//number//'.out',position='append')
        !        write(20,'(14E16.8)') 1.*istep,time/tref,ap(k)%x/lref,ap(k)%y/lref,ap(k)%z/lref, &
        !                             ap(k)%theta,ap(k)%phi,ap(k)%u/uref,ap(k)%v/uref,ap(k)%w/uref, &
        !                             ap(k)%omx*tref,ap(k)%omy*tref,ap(k)%omz*tref, &
        !                             ap(k)%omtheta*tref
        !        close(20)
        !      endif
        !    enddo
        call post1d(time,istep,nstep)  
        !   if (mod(istep,500).eq.0) then  ! (COBP)
      !   if (np > 0) then ! not needed at this point (initial check single particle FF in homogeneous) ! (SS)
      !      call autocorr(istep) ! Restore setup in Pedro's original files (SS)
      !      !   end if  ! (COBP)
      !   end if
     endif

     if (mod(istep,ioutchk).eq.0) then
        if (np > 0) then
           call loadpart(1,istep)           
        end if  
      !   call loadflds(1,istep)
     end if

     if (mod(istep,iout2d).eq.0) then
        call post3d(istep)
     end if     
     ! output routines


     ! Not needed at this point ! (MVR)
   !   if (mod(istep,iout2d).eq.0) then
   !      call post2d(istep)
   !      !call post3d(istep)  ! (COBP) in the next 6 lines
   !      !if(ihistory.eq.0) then
   !      !  count_history_max = nint(lz/bulk_v_sup/dt)
   !      !  ihistory = 1
   !      !  count_history = 0
   !      !  ihistory_acc = ihistory_acc + 1
   !      !endif
   !   endif
   !   !  if(ihistory.eq.1) then ! (COBP) in the next 5 lines
   !   !    count_history = count_history + 1
   !   !    call historypart(ihistory_acc,count_history)
   !   !    if(count_history.gt.count_history_max) ihistory = 0
   !   !  endif
   !   !
   !   if (mod(istep,ioutfld).eq.0) then
   !      !    call loadflds(1,istep) ! (COBP)
   !      !    call loadpart(1,istep) ! (COBP)
   !   endif

     ttemp = MPI_WTIME()-tsave
     call mpi_allreduce(ttemp,comp_time_avrg,1,mpi_real8,mpi_sum,comm_cart,error)
     ttemp = comp_time_avrg/product(dims)
     if (ttemp.gt.tmax*3600..or.istep.eq.nstep.or.time.ge.t_end) then
        !  if (mod(istep,ioutfld).eq.0) then ! (COBP) in the next 2 lines
        !    call loadpart(1,istep)
        !    call loadflds(1,istep)
        if (time.ge.t_end .and. myid .eq. 0) then
           write(6,*) 'Simulation stopped: time reached', t_end
        endif
        if (np > 0) then
           call loadpart(1,istep)
        end if
        call loadflds(1,istep)
        call mpi_barrier(comm_cart,error)
        goto 111
     endif
     !
     ! timming
     !
     comp_time = MPI_WTIME()-comp_time
     call mpi_allreduce(comp_time,comp_time_avrg,1,mpi_real8,mpi_sum,comm_cart,error)
     call mpi_allreduce(comp_time,comp_time_min,1,mpi_real8,mpi_min,comm_cart,error)
     call mpi_allreduce(comp_time,comp_time_max,1,mpi_real8,mpi_max,comm_cart,error)
     if(myid.eq.0) write(6,'(A,3E16.8)') 'Avrg, min & max elapsed time = ', &
          comp_time_avrg/product(dims),comp_time_min,comp_time_max

    
  enddo

111 continue
  !
  ! clean fftw used for autocorrelations
  !
  call clean_fft(plan_r2c_x,plan_c2r_x)
  call clean_fft(plan_r2c_y,plan_c2r_y)
  !
  if(myid.eq.0) write(6,*) '*** Fim ***'
  !
  call decomp_2d_finalize ! release resources used by decomp_2d
  call MPI_FINALIZE(error) ! finalize MPI
  !
  stop
end program interpt
