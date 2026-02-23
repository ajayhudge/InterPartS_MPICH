module mod_init
  use mod_common
  use decomp_2d
  use mod_param
  use mod_common_mpi
  implicit none
  private
  public init_ics, input_parameters, init_grid, derived_parameters
contains
  SUBROUTINE input_parameters
    integer:: ios, ierr
    logical:: lexist
    character(len=30) :: fname

    fname = "parameters.inp"

    inquire(file=fname, EXIST=lexist)
    if ( .not. lexist) then
      if (myid == 0) write(*, *) 'parameters NOT found, exiting ..'
      call decomp_2d_finalize
      call MPI_FINALIZE(error)
      STOP
    else
      
      if( myid == 0 ) then

        open(11,file=fname,status='old',iostat=ios)
        if( ios /= 0 ) then 
          write(*,*) "parameters not found, exiting .."
          call MPI_Abort(MPI_COMM_WORLD,MPI_ERR_OTHER,ierr)
          STOP
        else
          ! read the namelist namrun, defined in parameters.f90
          read(11,parameters)
          close(11)
        end if
      endif
    endif

   call mpi_bcast(lx, 1, MPI_DOUBLE_PRECISION, 0, comm_cart, ierr)
   call mpi_bcast(ly, 1, MPI_DOUBLE_PRECISION, 0, comm_cart, ierr)
   call mpi_bcast(lz, 1, MPI_DOUBLE_PRECISION, 0, comm_cart, ierr)
   call mpi_bcast(iniu, len(iniu), MPI_CHARACTER, 0, comm_cart, ierr)
   call mpi_bcast(inic, len(inic), MPI_CHARACTER, 0, comm_cart, ierr)
   call mpi_bcast(bc_uvw_top_type, len(bc_uvw_top_type), MPI_CHARACTER, 0, comm_cart, ierr)
   call mpi_bcast(bc_uvw_bot_type, len(bc_uvw_bot_type), MPI_CHARACTER, 0, comm_cart, ierr)
   call mpi_bcast(bc_c_top_type, len(bc_c_top_type), MPI_CHARACTER, 0, comm_cart, ierr)
   call mpi_bcast(bc_c_bot_type, len(bc_c_bot_type), MPI_CHARACTER, 0, comm_cart, ierr)
   call mpi_bcast(bc_c_top_val, 1, MPI_DOUBLE_PRECISION, 0, comm_cart, ierr)
   call mpi_bcast(bc_c_bot_val, 1, MPI_DOUBLE_PRECISION, 0, comm_cart, ierr)
   call mpi_bcast(particle_rotation, 1, MPI_LOGICAL, 0, comm_cart, ierr)
   ! call mpi_bcast(isfreeslip, 1, MPI_LOGICAL, 0, comm_cart, ierr)
   call mpi_bcast(visc, 1, MPI_DOUBLE_PRECISION, 0, comm_cart, ierr)
   call mpi_bcast(Pran, 1, MPI_DOUBLE_PRECISION, 0, comm_cart, ierr)
   call mpi_bcast(beta, 1, MPI_DOUBLE_PRECISION, 0, comm_cart, ierr)
   call mpi_bcast(radius, 1, MPI_DOUBLE_PRECISION, 0, comm_cart, ierr)
   call mpi_bcast(ratiorho, 1, MPI_DOUBLE_PRECISION, 0, comm_cart, ierr)
   call mpi_bcast(ioutchk, 1, MPI_INTEGER, 0, comm_cart, ierr)
   call mpi_bcast(iout1d, 1, MPI_INTEGER, 0, comm_cart, ierr)
   call mpi_bcast(iout2d, 1, MPI_INTEGER, 0, comm_cart, ierr)
   call mpi_bcast(ioutfld, 1, MPI_INTEGER, 0, comm_cart, ierr)
   call mpi_bcast(runmode, 1, MPI_INTEGER, 0, comm_cart, ierr)

     ! let last processor report the contents of the namelist
    IF( myid == nproc - 1 ) THEN
      write(6,*) 'info from process', myid
      write(6,parameters)
    ENDIF

    call derived_parameters
  END SUBROUTINE input_parameters

  ! ------------------------------------------------------------------------- !
  SUBROUTINE derived_parameters
    ! This subroutine calculates derived parameters based on the input parameters.
    ! It is called after input_parameters to ensure all parameters are set.

    ! Grid properties
    dxi = itot/lx; dyi = jtot/ly; dzi = ktot/lz
    dx = 1./dxi; dy = 1./dyi; dz = 1./dzi
    
    
    ! boundaries of the local subdomain (SS)
    boundleftmyid  = coords(1)*lx/(1.*dims(1)) ! left  boundary
    boundfrontmyid = coords(2)*ly/(1.*dims(2)) ! front boundary



    ! amount of data to be sent from mstr to slve: see common file
    send_real = 24+10*nqmax
    send_int = 1+nqmax 

    ! Thermal diffusivity
    kappa = visc/Pran

    ! radius = 0.5 ! set particle radius as a parameter (SS)
    offset = ( sqrt(3.*(1.5**2)) )/dxi + 0.01/dxi
    volp = (4./3.)*pi*radius**3.
    mominert = (2./5.)*volp*radius**2.
    lref = 1. ! lz
    uref = 1. ! bulk_v_sup
    tref = 1. !lref/uref

    ! These checks are only relevant if particles are present
    if ( np > 0) then
       if ( (radius+offset) .gt. lx/2./dims(1)) then
          if (myid.eq.0) then
             write(6,*) 'Diameter spheres larger than x-dimension of processes'
             write(6,*) 'Change the value of dims(1) in "param.f90".'
             write(6,*) 'Program aborted...'
          endif
          call mpi_finalize(error)
          stop
       endif
       if ( (radius+offset) .gt. ly/2./dims(2)) then
          if (myid.eq.0) then
             write(6,*) 'Diameter spheres larger than y-dimension of processes'
             write(6,*) 'Change the value of dims(2) in "param.f90".'
             write(6,*) 'Program aborted...'
          endif
          call mpi_finalize(error)
          stop
       endif
    endif

    ! Particle properties 
    retrac = 0.3/dxi
  
    nl = nint((pi/3.)*(12.*(((radius-retrac)*dxi)**2)+1.)) 
    NL2=nint((pi/3.)*(12.*(((radius-1.5/dxi)*dxi)**2.) + 1. )) ! nr lfp's of second shell
    NL3=nint((pi/3.)*(12.*(((radius-2.5/dxi)*dxi)**2.) + 1. )) ! nr lfp's of third shell
    NL4=nint((pi/3.)*(12.*(((radius-3.5/dxi)*dxi)**2.) + 1. )) ! nr lfp's of fourth shell
    NLtot=NL+NL2+NL3+NL4 ! nr lfp's of 1st-4th shell
    
    solidity = (1.*np)*(4./3.)*pi*(radius**3)/(lx*ly*lz)
  
    ! sphere/sphere
    colthr_pp = 0.*radius
    meffn_ss = ratiorho*volp/2.
    mefft_ss = 2./7.*meffn_ss
    kn_ss = (pi**2. + abs(log(en))**2.)*meffn_ss/((Nstretch*dt_estim)**2.)
    kt_ss = (pi**2. + abs(log(et))**2.)*mefft_ss/((Nstretch*dt_estim)**2.)
    etan_ss = -2.*(log(en))*meffn_ss/(Nstretch*dt_estim)
    etat_ss = -2.*(log(et))*mefft_ss/(Nstretch*dt_estim)
    muc_ss = muc
    psi_crit_ss = 7./2.*(1.+en)/(1.+et)*muc_ss
  
    ! sphere/wall
    colthr_pw = 0.001*radius
    meffn_sw = ratiorho*volp
    mefft_sw = 2./7.*meffn_sw
    kn_sw = (pi**2. + abs(log(en))**2.)*meffn_sw/((Nstretch*dt_estim)**2.)
    kt_sw = (pi**2. + abs(log(et))**2.)*mefft_sw/((Nstretch*dt_estim)**2.)
    etan_sw = -2.*(log(en))*meffn_sw/(Nstretch*dt_estim)
    etat_sw = -2.*(log(et))*mefft_sw/(Nstretch*dt_estim)
    muc_sw = muc
    psi_crit_sw = 7./2.*(1.+en)/(1.+et)*muc_sw
  
    !
    ! set parameters for the lubrication model
    !
    coeff_f = 6.*pi*visc*radius
    coeff_t = 8.*pi*visc*radius**2.

    ! Correction with Stokes amplification factor is added for values of epsilon smaller than
    ! eps_ini_pp
    a11_ini_pp = -1./4.*eps_ini_pp**(-1.)+9./40.*log(eps_ini_pp)+3./112.*eps_ini_pp*log(eps_ini_pp)
    a11a_ini_pp = a11_ini_pp-.995
    a11b_ini_pp = -a11_ini_pp+.350
    a22_ini_pp = 1./6.*log(eps_ini_pp)
    a22a_ini_pp = a22_ini_pp-.998
    a22b_ini_pp = -a22_ini_pp+.274
    a33a_ini_pp = a22a_ini_pp
    a33b_ini_pp = a22b_ini_pp
    b23_ini_pp = -1./6.*log(eps_ini_pp)-1./12.*eps_ini_pp*log(eps_ini_pp)
    b23a_ini_pp = b23_ini_pp-.159
    b23b_ini_pp = -b23_ini_pp+.001
    b32a_ini_pp = -b23a_ini_pp
    b32b_ini_pp = -b23b_ini_pp
    c23a_ini_pp = b32a_ini_pp
    c23b_ini_pp = b32b_ini_pp
    c32a_ini_pp = b23a_ini_pp
    c32b_ini_pp = b23b_ini_pp
    d11_ini_pp = 1./8.*eps_ini_pp*log(eps_ini_pp)
    d11a_ini_pp = 1./8.*eps_ini_pp*log(eps_ini_pp)
    d11b_ini_pp = -1./8.*eps_ini_pp*log(eps_ini_pp)
    d22a_ini_pp = 1./5.*log(eps_ini_pp)+47./250.*eps_ini_pp*log(eps_ini_pp)-.703
    d22b_ini_pp = -1./20.*log(eps_ini_pp)+31./500.*eps_ini_pp*log(eps_ini_pp)-.027
    d33a_ini_pp = 1./5.*log(eps_ini_pp)+47./250.*eps_ini_pp*log(eps_ini_pp)-.703
    d33b_ini_pp = -1./20.*log(eps_ini_pp)+31./500.*eps_ini_pp*log(eps_ini_pp)-.027
  
    ! Stokes amplification factor saturated for values of epsilon smaller than eps_sat_pp
    a11_sat_pp = -1./4.*eps_sat_pp**(-1.)+9./40.*log(eps_sat_pp)+3./112.*eps_sat_pp*log(eps_sat_pp)
    a11a_sat_pp = a11_sat_pp-.995
    a11b_sat_pp = -a11_sat_pp+.350
    a22_sat_pp = 1./6.*log(eps_sat_pp)
    a22a_sat_pp = a22_sat_pp-.998
    a22b_sat_pp = -a22_sat_pp+.274
    a33a_sat_pp = a22a_sat_pp
    a33b_sat_pp = a22b_sat_pp
    b23_sat_pp = -1./6.*log(eps_sat_pp)-1./12.*eps_sat_pp*log(eps_sat_pp)
    b23a_sat_pp = b23_sat_pp-.159
    b23b_sat_pp = -b23_sat_pp+.001
    b32a_sat_pp = -b23a_sat_pp
    b32b_sat_pp = -b23b_sat_pp
    c23a_sat_pp = b32a_sat_pp
    c23b_sat_pp = b32b_sat_pp
    c32a_sat_pp = b23a_sat_pp
    c32b_sat_pp = b23b_sat_pp
    d11_sat_pp = 1./8.*eps_sat_pp*log(eps_sat_pp)
    d11a_sat_pp = 1./8.*eps_sat_pp*log(eps_sat_pp)
    d11b_sat_pp = -1./8.*eps_sat_pp*log(eps_sat_pp)
    d22a_sat_pp = 1./5.*log(eps_sat_pp)+47./250.*eps_sat_pp*log(eps_sat_pp)-.703
    d22b_sat_pp = -1./20.*log(eps_sat_pp)+31./500.*eps_sat_pp*log(eps_sat_pp)-.027
    d33a_sat_pp = 1./5.*log(eps_sat_pp)+47./250.*eps_sat_pp*log(eps_sat_pp)-.703
    d33b_sat_pp = -1./20.*log(eps_sat_pp)+31./500.*eps_sat_pp*log(eps_sat_pp)-.027

    ! Correction with Stokes amplification factor is added for values of epsilon smaller than
    ! eps_ini_pw
    a11_ini_pw = -1./eps_ini_pw+1./5.*log(eps_ini_pw)+1./21.*eps_ini_pw*log(eps_ini_pw)-.9713
    a22_ini_pw = 8./15.*log(eps_ini_pw)+64./375.*eps_ini_pw*log(eps_ini_pw)-.952
    a33_ini_pw = a22_ini_pw
    b23_ini_pw = -2./15.*log(eps_ini_pw)-86./375.*eps_ini_pw*log(eps_ini_pw)-.257
    b32_ini_pw = -b23_ini_pw
    c23_ini_pw = b32_ini_pw
    c32_ini_pw = b23_ini_pw
    d11_ini_pw = 1./2.*eps_ini_pw*log(eps_ini_pw)-1.202
    d22_ini_pw = 2./5.*log(eps_ini_pw)+66./125.*eps_ini_pw*log(eps_ini_pw)-.371
    d33_ini_pw = 2./5.*log(eps_ini_pw)+66./125.*eps_ini_pw*log(eps_ini_pw)-.371
  
    ! Stokes amplification factor saturated for values of epsilon smaller than eps_sat_pw
    a11_sat_pw = -1./eps_sat_pw+1./5.*log(eps_sat_pw)+1./21.*eps_sat_pw*log(eps_sat_pw)-.9713
    a22_sat_pw = 8./15.*log(eps_sat_pw)+64./375.*eps_sat_pw*log(eps_sat_pw)-.952
    a33_sat_pw = a22_sat_pw
    b23_sat_pw = -2./15.*log(eps_sat_pw)-86./375.*eps_sat_pw*log(eps_sat_pw)-.257
    b32_sat_pw = -b23_sat_pw
    c23_sat_pw = b32_sat_pw
    c32_sat_pw = b23_sat_pw
    d11_sat_pw = 1./2.*eps_sat_pw*log(eps_sat_pw)-1.202
    d22_sat_pw = 2./5.*log(eps_sat_pw)+66./125.*eps_sat_pw*log(eps_sat_pw)-.371
    d33_sat_pw = 2./5.*log(eps_sat_pw)+66./125.*eps_sat_pw*log(eps_sat_pw)-.371

  END SUBROUTINE derived_parameters

  ! ------------------------------------------------------------------------- !
  !
  SUBROUTINE init_grid
    integer :: i,j,k

    DO i = 0, i1
      xb(i) = (i + coords(1)*imax)*dx ! x coordinate for face centred variable
      xc(i) = (i - 0.5 + coords(1)*imax)*dx  ! x coordinate for cell centred variable
    ENDDO

    DO j = 0, j1
      yb(j) = (j + coords(2)*jmax)*dy
      yc(j) = (j - 0.5 + coords(2)*jmax)*dy
    ENDDO

    DO k = 0, k1
      zb(k) = k*dz
      zc(k) = (k-0.5)*dz
    ENDDO    
    
    ! call showgrid
  END SUBROUTINE

  ! ------------------------------------------------------------------------- !  

  subroutine init_ics
    implicit none
    integer :: i,j,k,iglob,jglob,kglob
    integer, allocatable :: seed(:)
    real :: coorz,coory,coorx,rn1,rn2,rn3
    real, dimension(kmax) :: v
    real umean,vmean,wmean
    real umean_all,vmean_all,wmean_all
    real, parameter :: bulk_v_sup = 1.0
    real :: Re

    ! initial conditions
    !
    select case(iniu)
       !
    case('zer') ! velocity field = 0 forall i,j,k with small noise
       ! Initialize random seed differently for each MPI rank
       call random_seed()
       ! Use rank-dependent seed to ensure different random numbers on each process
       call random_seed(size=k)  ! Get size of seed array
       allocate(seed(k))
       call random_seed(get=seed)
       seed = seed + myid * 1000  ! Make seed unique for each rank
       call random_seed(put=seed)
       deallocate(seed)
       
       do k=0,k1
          do j=0,j1
             do i=0,i1
                call random_number(rn1)
                call random_number(rn2)
                call random_number(rn3)
                ! Add noise with magnitude in range [1.5e-4, 3.3e-4] and random sign
                unew(i,j,k) = (1.5e-3 + rn1 * 1.8e-3) * sign(1.0, rn1 - 0.5)
                vnew(i,j,k) = (1.5e-3 + rn2 * 1.8e-3) * sign(1.0, rn2 - 0.5)
                wnew(i,j,k) = (1.5e-3 + rn3 * 1.8e-3) * sign(1.0, rn3 - 0.5)
                ! dudt(i,j,k) = 0.
                ! dvdt(i,j,k) = 0.
                ! dwdt(i,j,k) = 0.
                pnew(i,j,k) = 0.
             enddo
          enddo
       enddo
       !
    case('pur') ! pure zero velocity field
       do k=0,k1
          do j=0,j1
             do i=0,i1
                unew(i,j,k) = 0.0
                vnew(i,j,k) = 0.0
                wnew(i,j,k) = 0.0
                ! dudt(i,j,k) = 0.
                ! dvdt(i,j,k) = 0.
                ! dwdt(i,j,k) = 0.
                pnew(i,j,k) = 0.
             enddo
          enddo
       enddo
       !       
    case('cou') ! Couette flow
       do k=1,kmax
          coorz = (k-0.5)/dzi/lz ! normalised with channel height
          do j=1,jmax
             do i=1,imax
                unew(i,j,k) = 0.
                vnew(i,j,k) = .5*bulk_v_sup*(1.-2.*coorz)
                wnew(i,j,k) = 0.
                ! dudt(i,j,k) = 0.
                ! dvdt(i,j,k) = 0.
                ! dwdt(i,j,k) = 0.
                pnew(i,j,k) = 0.
             enddo
          enddo
       enddo
       !
    case('log') ! Logarithmic profile
      !  if(isfreeslip) then
      !     do k=1,kmax
      !        v(k) = 2.5*log( (1./visc)*(k-0.5)/(2.*kmax) ) + 5.5
      !        if ((k-0.5)/kmax/visc .le. 11.6) v(k)=(k-0.5)/kmax/visc
      !     enddo
      !  else
      ! Shuang: The reason for commenting out the above lines is that I have replaced 'isfreeslip' 
      ! with variables like BC_NOSLIP in param.f90. However, I am not sure about the exact 
      ! role of 'isfreeslip' here, so I commented it out for now.
      do k = 1, kmax/2
         v(k) = 2.5*log((1./visc)*(k-0.5)/(1.*kmax)) + 5.5
         if ((k-0.5)/kmax/visc .le. 11.6) v(k) = (k-0.5)/kmax/visc
      enddo
      do k = 1, kmax/2
         i = kmax + 1 - k
         v(kmax+1-k) = 2.5*log((1./visc)*(k-0.5)/(1.*kmax)) + 5.5
         if ((k-0.5)/kmax/visc .le. 11.6) v(kmax+1-k) = (k-0.5)/kmax/visc
      enddo
      !  endif
       !
       ! Below, random numbers are generated such that the initial field
       ! does not depend in the number of mpi tasks.
       !
       umean=0.
       vmean=0.
       wmean=0.
       !  call random_seed( put = (/16112006/) )
       !  do kglob=1,ktot
       !    do jglob=1,jtot
       !      do iglob=1,itot
       !        call random_number(rn1)
       !        call random_number(rn2)
       !        call random_number(rn3)
       !        i = iglob-imax*coords(1)
       !        j = jglob-jmax*coords(2)
       !        k = kglob
       !        if( &
       !            (1.le.i.and.i.le.imax) .and. &
       !            (1.le.j.and.j.le.jmax) &
       !          ) then
       !          unew(i,j,k)=15.*(rn1-0.5)
       !          umean=umean+unew(i,j,k)
       !          vnew(i,j,k)=v(k)*(1.+5.*(rn2-0.5))
       !          vmean=vmean+vnew(i,j,k)
       !          wnew(i,j,k)=15.*(rn3-0.5)
       !          wmean=wmean+wnew(i,j,k)
       !          dudt(i,j,k)=0.
       !          dvdt(i,j,k)=0.
       !          dwdt(i,j,k)=0.
       !          pnew(i,j,k)=0.
       !        endif
       !      enddo
       !    enddo
       !  enddo
       do k=1,kmax
          do j=1,jmax
             do i=1,imax
                call random_number(rn1)
                unew(i,j,k)=15.*(rn1-0.5)
                umean=umean+unew(i,j,k)
                call random_number(rn1)
                vnew(i,j,k)=v(k)*(1.+5.*(rn1-0.5))
                vmean=vmean+vnew(i,j,k)
                call random_number(rn1)
                wnew(i,j,k)=15.*(rn1-0.5)
                wmean=wmean+wnew(i,j,k)
                ! dudt(i,j,k)=0.
                ! dvdt(i,j,k)=0.
                ! dwdt(i,j,k)=0.
                pnew(i,j,k)=0.
             enddo
          enddo
       enddo
       call mpi_allreduce(umean,umean_all,1,mpi_real8,mpi_sum,comm_cart,error)
       call mpi_allreduce(vmean,vmean_all,1,mpi_real8,mpi_sum,comm_cart,error)
       call mpi_allreduce(wmean,wmean_all,1,mpi_real8,mpi_sum,comm_cart,error)
       umean_all = umean_all/(1.*itot*jtot*kmax)
       vmean_all = vmean_all/(1.*itot*jtot*kmax)
       wmean_all = wmean_all/(1.*itot*jtot*kmax)
       
       do k=1,kmax
          do j=1,jmax
             do i=1,imax
                unew(i,j,k) = unew(i,j,k)-umean_all !average = 0
                vnew(i,j,k) = vnew(i,j,k)/vmean_all !average = 1
                wnew(i,j,k) = wnew(i,j,k)-wmean_all !average = 0
             enddo
          enddo
       enddo
       !
    case('poi') ! Poiseuille
       do k = 0,k1
          coorz = (k-0.5)*dz/lz ! normalised with channel height
          do j = 0,j1
             do i = 0,i1
                unew(i,j,k) = 0.
                vnew(i,j,k) = 6.*bulk_v_sup*coorz*(1.-coorz)
                wnew(i,j,k) = 0.
                ! dudt(i,j,k) = 0.
                ! dvdt(i,j,k) = 0.
                ! dwdt(i,j,k) = 0.
                pnew(i,j,k) = 0.
             enddo
          enddo
       enddo
 
    case('TGV') ! Taylor-Green Vortex in x and z directions, decaying in time
      ! Converge test, need to check whether coords(1) affect the result (SS)
       do k = 1,kmax
          do j = 1,jmax
             do i = 1,imax
                unew(i,j,k) =  sin(2*pi*xb(i)) * cos(2*pi*zc(k))
                vnew(i,j,k) = 0.0
                wnew(i,j,k) = -cos(2*pi*xc(i)) * sin(2*pi*zb(k))
             end do
          end do
       end do

    end select


    select case(inic)
    !  
    case('RBC') ! Rayleigh-Bénard Convection
       do k = 1,kmax
          do j = 1,jmax
             do i = 1,imax
                cnew(i,j,k) = bc_c_top_val + (bc_c_bot_val - bc_c_top_val)* (1 - zc(k)/lz) + &
                epsilon*sin(pi*zc(k)/lz)*sin(pi*xc(i)/lx)*sin(pi*yc(j)/ly)
                ! linear distribution + small perturbation in the domain
             end do
          end do
       end do
    !  
    case('ZERC') ! Zero concentration field
       do k = 1,kmax
          do j = 1,jmax
             do i = 1,imax
                cnew(i,j,k) = 0.0
             end do
          end do
       end do
    !

    !  
    case('FFS') ! Free Falling Sphere with stratified fluids controled by temperature
       do k = 1,kmax
          do j = 1,jmax
             do i = 1,imax
                cnew(i,j,k) = (bc_c_top_val + bc_c_bot_val)/2. + (bc_c_top_val - bc_c_bot_val)/2*tanh((zc(k) - 0.37*lz)*4/h_inter) 
                ! linear distribution + small perturbation in the domain
             end do
          end do
       end do
    !  

    end select
    !
    return
  end subroutine init_ics
end module mod_init

          
