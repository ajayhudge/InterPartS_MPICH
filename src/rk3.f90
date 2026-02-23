module mod_rk
   use mod_mom
   use mod_scal
   use mod_common
   use mod_common_mpi
   use mod_test_sources
   use mod_param, only: runmode
  implicit none
  private
   public rk1,rk2,rk3,rk_stage
contains
  !
  subroutine rk1(durk1,dvrk1,dwrk1,dcrk1)
    !
    ! First step of a low-storage 3rd-order Runge-Kutta scheme 
    ! for time integration of the momentum equations.
    ! Out : rhs of mom. eqs at RK1 level
    !
    implicit none
    integer i,j,k
    real, intent(inout) :: durk1(0:,0:,0:)
    real, intent(inout) :: dvrk1(0:,0:,0:)
    real, intent(inout) :: dwrk1(0:,0:,0:)
    real, intent(inout) :: dcrk1(0:,0:,0:)
    real :: dwrk1_old(0:i1,0:j1,0:k1)
    real :: wallshear_all
    real :: factor
    real :: dVeul
    !
    factor = dt*(32./60.)
    !
    call momxp(durk1,pnew)
    call momyp(dvrk1,pnew)
    call momzp(dwrk1,pnew)
    !$omp parallel default(none)    & 
    !$omp&shared(factor)            &
    !$omp&shared(durk1,dvrk1,dwrk1,dcrk1) &
    !$omp&shared(unew ,vnew ,wnew,cnew ) &
    !$omp&shared(dudt ,dvdt ,dwdt,dcdt) &
    !$omp&private(i,j,k)
    !$omp do
    do k=1,kmax
       do j=1,jmax
          do i=1,imax
             dudt(i,j,k) = unew(i,j,k) + factor*durk1(i,j,k)
             dvdt(i,j,k) = vnew(i,j,k) + factor*dvrk1(i,j,k)
             dwdt(i,j,k) = wnew(i,j,k) + factor*dwrk1(i,j,k)
             dcdt(i,j,k) = cnew(i,j,k) ! + factor*dcrk1(i,j,k) 
          enddo
       enddo
    enddo
    !$omp end parallel
    call momxad(durk1,unew,vnew,wnew)
    call momyad(dvrk1,unew,vnew,wnew)
    call momzad(dwrk1,unew,vnew,wnew)
    dwrk1_old(:,:,:) = dwrk1(:,:,:) ! (SS)
    call momzbuoy(dwrk1,cnew)
    dwrk1(:,:,:) = dwrk1(:,:,:) + dwrk1_old(:,:,:) ! (SS)
    call scalad(dcrk1,unew,vnew,wnew,cnew)


    if (runmode == RUNMODE_TEST_SCALAR_SCHEME) then
       call tests_scalar_scheme_src(dcrk1,dwrk1,time)
    end if

    !$omp parallel default(none)    &
    !$omp&shared(factor)            &
    !$omp&shared(durk1,dvrk1,dwrk1,dcrk1) &
    !$omp&shared(dudt ,dvdt ,dwdt,dcdt) &
    !$omp&private(i,j,k)
    !$omp do
    do k=1,kmax
       do j=1,jmax
          do i=1,imax
             dudt(i,j,k) = dudt(i,j,k) + factor*durk1(i,j,k)
             dvdt(i,j,k) = dvdt(i,j,k) + factor*dvrk1(i,j,k)
             dwdt(i,j,k) = dwdt(i,j,k) + factor*dwrk1(i,j,k)
             dcdt(i,j,k) = dcdt(i,j,k) + factor*dcrk1(i,j,k)
          enddo
       enddo
    enddo
    ! cnew(:,:,:) = dcdt(:,:,:) !(SS)
    !$omp end parallel
    !
    ! computation averaged wall shear stress force
    !
    dVeul = (1./dxi)*(1./dyi)*(1./dzi)
    wallshearnew = 0.
   !!$omp parallel default(shared) &
   !!$omp&private(i,j,k) &
   !!$omp&reduction(+:wallshearnew)
   !!$omp do
    do j=1,jmax
      do i=1,imax
        !no-slip condition at both boundaries
        wallshearnew = wallshearnew + visc*(vnew(i,j,kmax+1)-vnew(i,j,kmax))*dzi*(1./dxi)*(1./dyi) &
                                    - visc*(vnew(i,j,1)-vnew(i,j,0))*dzi*(1./dxi)*(1./dyi)
      enddo
    enddo
   !!$omp end parallel
    call mpi_allreduce(wallshearnew,wallshear_all,1,mpi_real8,mpi_sum,comm_cart,error)
    wallshearnew = (factor/dt)*wallshear_all/(dVeul*itot*jtot*kmax)
    wallshearold = wallshear_all/(dVeul*itot*jtot*kmax)
    !!
    return
  end subroutine rk1
  !
  subroutine rk2(durk1,dvrk1,dwrk1,dcrk1)
    !
    ! Second step of a low-storage 3rd-order Runge-Kutta scheme 
    ! for time integration of the momentum equations.
    ! In  : rhs of mom. eqs at RK1 level
    ! Out : rhs of mom. eqs at RK2 level
    !
    implicit none
    integer i,j,k
    real, intent(inout) :: durk1(0:,0:,0:)
    real, intent(inout) :: dvrk1(0:,0:,0:)
    real, intent(inout) :: dwrk1(0:,0:,0:)
    real, intent(inout) :: dcrk1(0:,0:,0:)
    real :: dnewu(0:i1,0:j1,0:k1) ! dummy array
    real :: dnewv(0:i1,0:j1,0:k1) ! dummy array
    real :: dneww(0:i1,0:j1,0:k1) ! dummy array
    real :: dnewc(0:i1,0:j1,0:k1) ! dummy array
    real :: dneww_old(0:i1,0:j1,0:k1)
    real :: wallshear_all
    real :: factor1,factor2,factor3
    real :: dVeul
    character(len=64) :: fname3,fname4
    !
    factor1 = dt*( 8./60.)
    factor2 = dt*(25./60.)
    factor3 = dt*(17./60.)
    !
    call momxp(dnewu,pnew)
    call momyp(dnewv,pnew)
    call momzp(dneww,pnew)
    !$omp parallel default(none)    &
    !$omp&shared(factor1)           &
    !$omp&shared(dnewu,dnewv,dneww,dnewc) &
    !$omp&shared(unew ,vnew ,wnew,cnew ) &
    !$omp&shared(dudt ,dvdt ,dwdt,dcdt) &
    !$omp&private(i,j,k)
    !$omp do
    do k=1,kmax
       do j=1,jmax
          do i=1,imax
             dudt(i,j,k) = unew(i,j,k) + factor1*dnewu(i,j,k)
             dvdt(i,j,k) = vnew(i,j,k) + factor1*dnewv(i,j,k)
             dwdt(i,j,k) = wnew(i,j,k) + factor1*dneww(i,j,k)
             dcdt(i,j,k) = cnew(i,j,k) ! + factor1*dnewc(i,j,k) 
          enddo
       enddo
    enddo
    !$omp end parallel
    call momxad(dnewu,unew,vnew,wnew)
    call momyad(dnewv,unew,vnew,wnew)
    call momzad(dneww,unew,vnew,wnew)
    dneww_old(:,:,:) = dneww(:,:,:) ! (SS)
    call momzbuoy(dneww,cnew)
    dneww(:,:,:) = dneww(:,:,:) + dneww_old(:,:,:) ! (SS)
    call scalad(dnewc,unew,vnew,wnew,cnew) ! scalar



    if (runmode == RUNMODE_TEST_SCALAR_SCHEME) then
       call tests_scalar_scheme_src(dnewc,dneww,time)
    end if



    !$omp parallel default(none)    &
    !$omp&shared(factor2,factor3)   &
    !$omp&shared(durk1,dvrk1,dwrk1,dcrk1) &
    !$omp&shared(dnewu,dnewv,dneww,dnewc) &
    !$omp&shared(dudt ,dvdt ,dwdt ,dcdt) &
    !$omp&private(i,j,k)
    !$omp do
    do k=1,kmax
       do j=1,jmax
          do i=1,imax
             dudt(i,j,k) = dudt(i,j,k) + & 
                          ( factor2*dnewu(i,j,k) - factor3*durk1(i,j,k) )
             dvdt(i,j,k) = dvdt(i,j,k) + &
                          ( factor2*dnewv(i,j,k) - factor3*dvrk1(i,j,k) )
             dwdt(i,j,k) = dwdt(i,j,k) + &
                          ( factor2*dneww(i,j,k) - factor3*dwrk1(i,j,k) )
             dcdt(i,j,k) = dcdt(i,j,k) + &
                          ( factor2*dnewc(i,j,k) - factor3*dcrk1(i,j,k) )
             durk1(i,j,k) = dnewu(i,j,k)
             dvrk1(i,j,k) = dnewv(i,j,k)
             dwrk1(i,j,k) = dneww(i,j,k)
             dcrk1(i,j,k) = dnewc(i,j,k)
          enddo
       enddo
    enddo
    ! cnew(:,:,:) = dcdt(:,:,:) !(SS)
    !$omp end parallel
    !
    ! computation averaged wall shear stress force
    ! 
    dVeul = (1./dxi)*(1./dyi)*(1./dzi)
    wallshearnew = 0. !initialisation
    !!$omp parallel default(shared) &
    !!$omp&private(i,j,k) &
    !!$omp&reduction(+:wallshearnew)
    !!$omp do
    do j=1,jmax
      do i=1,imax
        !no-slip condition at both boundaries
        wallshearnew = wallshearnew + visc*(vnew(i,j,kmax+1)-vnew(i,j,kmax))*dzi*(1./dxi)*(1./dyi) &
                                    - visc*(vnew(i,j,1)-vnew(i,j,0))*dzi*(1./dxi)*(1./dyi)
      enddo
    enddo
    !!$omp end parallel
    call mpi_allreduce(wallshearnew,wallshear_all,1,mpi_real8,mpi_sum,comm_cart,error)
    wallshearnew = (factor2/dt)*wallshear_all/(dVeul*itot*jtot*kmax) - &
                   (factor3/dt)*wallshearold
    wallshearold = wallshear_all/(dVeul*itot*jtot*kmax)
    !
    return
  end subroutine rk2
  !
  subroutine rk3(durk2,dvrk2,dwrk2,dcrk2)
    !
    ! Third step of a low-storage 3rd-order Runge-Kutta scheme 
    ! for time integration of the momentum equations.
    ! In  : rhs of mom. eqs at RK2 level
    !
    implicit none
    integer i,j,k
    real, intent(inout) :: durk2(0:,0:,0:)
    real, intent(inout) :: dvrk2(0:,0:,0:)
    real, intent(inout) :: dwrk2(0:,0:,0:)
    real, intent(inout) :: dcrk2(0:,0:,0:)
    real :: dnewu(0:i1,0:j1,0:k1) ! dummy array
    real :: dnewv(0:i1,0:j1,0:k1) ! dummy array
    real :: dneww(0:i1,0:j1,0:k1) ! dummy array
    real :: dnewc(0:i1,0:j1,0:k1) ! dummy array
    real :: dneww_old(0:i1,0:j1,0:k1)
    real :: wallshear_all
    real :: factor1,factor2,factor3
    real :: dVeul
    character(len=64) :: fname5,fname6
    !
    factor1 = dt*(20./60.)
    factor2 = dt*(45./60.)
    factor3 = dt*(25./60.)
    !
    call momxp(dnewu,pnew)
    call momyp(dnewv,pnew)
    call momzp(dneww,pnew)
    !$omp parallel default(none)    &
    !$omp&shared(factor1)           &
    !$omp&shared(dnewu,dnewv,dneww,dnewc) &
    !$omp&shared(unew ,vnew ,wnew ,cnew ) &
    !$omp&shared(dudt ,dvdt ,dwdt ,dcdt) &
    !$omp&private(i,j,k)
    !$omp do
    do k=1,kmax
       do j=1,jmax
          do i=1,imax
             dudt(i,j,k) = unew(i,j,k) + factor1*dnewu(i,j,k)
             dvdt(i,j,k) = vnew(i,j,k) + factor1*dnewv(i,j,k)
             dwdt(i,j,k) = wnew(i,j,k) + factor1*dneww(i,j,k)
             dcdt(i,j,k) = cnew(i,j,k) ! + factor1*dnewc(i,j,k)
          enddo
       enddo
    enddo
    !$omp end parallel
    call momxad(dnewu,unew,vnew,wnew)
    call momyad(dnewv,unew,vnew,wnew)
    call momzad(dneww,unew,vnew,wnew)
    dneww_old(:,:,:) = dneww(:,:,:) ! (SS)
    call momzbuoy(dneww,cnew)
    dneww(:,:,:) = dneww(:,:,:) + dneww_old(:,:,:) ! (SS)
    call scalad(dnewc,unew,vnew,wnew,cnew) ! scalar


    if (runmode == RUNMODE_TEST_SCALAR_SCHEME) then
       call tests_scalar_scheme_src(dnewc,dneww,time)
    end if

    !$omp parallel default(none)    &
    !$omp&shared(factor2,factor3)   &
    !$omp&shared(durk2,dvrk2,dwrk2,dcrk2) &
    !$omp&shared(dnewu,dnewv,dneww,dnewc) &
    !$omp&shared(dudt ,dvdt ,dwdt ,dcdt) &
    !$omp&private(i,j,k)
    !$omp do
    do k=1,kmax
       do j=1,jmax
          do i=1,imax
             dudt(i,j,k) = dudt(i,j,k) + & 
                           ( factor2*dnewu(i,j,k) - factor3*durk2(i,j,k) )
             dvdt(i,j,k) = dvdt(i,j,k) + &
                           ( factor2*dnewv(i,j,k) - factor3*dvrk2(i,j,k) )
             dwdt(i,j,k) = dwdt(i,j,k) + &
                           ( factor2*dneww(i,j,k) - factor3*dwrk2(i,j,k) )
             dcdt(i,j,k) = dcdt(i,j,k) + &
                           ( factor2*dnewc(i,j,k) - factor3*dcrk2(i,j,k) )
          enddo
       enddo
    enddo
    ! cnew(:,:,:) = dcdt(:,:,:) !(SS)
    !$omp end parallel
    !
    ! computation averaged wall shear stress force
    !
    dVeul = (1./dxi)*(1./dyi)*(1./dzi)
    wallshearnew = 0. !initialisation
   !!$omp parallel default(shared) &
   !!$omp&private(i,j,k) &
   !!$omp&reduction(+:wallshearnew)
   !!$omp do
    do j=1,jmax
      do i=1,imax
        !no-slip condition at both boundaries
        wallshearnew = wallshearnew + visc*(vnew(i,j,kmax+1)-vnew(i,j,kmax))*dzi*(1./dxi)*(1./dyi) &
                                    - visc*(vnew(i,j,1)-vnew(i,j,0))*dzi*(1./dxi)*(1./dyi)
      enddo
    enddo
   !!$omp end parallel
    call mpi_allreduce(wallshearnew,wallshear_all,1,mpi_real8,mpi_sum,comm_cart,error)
    wallshearnew = (factor2/dt)*wallshear_all/(dVeul*itot*jtot*kmax) - &
                   (factor3/dt)*wallshearold
    !              
    return
  end subroutine rk3
  !
   ! -----------------------------------------------------------------------
   ! Unified RK stage driver (not wired in yet)
   ! pass = 0 -> rk1, pass = 1 -> rk2, pass = 2 -> rk3
   subroutine rk_stage(pass, durk, dvrk, dwrk, dcrk)
      use mod_mom
      use mod_scal
      use mod_common
      use mod_common_mpi
      use mod_param
      implicit none
      integer, intent(in) :: pass
      real, dimension(:,:,:), intent(inout) :: durk, dvrk, dwrk, dcrk
      ! locals
      integer :: i,j,k
      real :: dnewu(0:i1,0:j1,0:k1)
      real :: dnewv(0:i1,0:j1,0:k1)
      real :: dneww(0:i1,0:j1,0:k1)
      real :: dnewc(0:i1,0:j1,0:k1)
      real :: dwrk_old(0:i1,0:j1,0:k1)
      real :: dneww_old(0:i1,0:j1,0:k1)
      real :: f1, f2, f3
      logical :: use_local, copy_back
   real :: wallshear_all, dVeul_loc, avgShear

      select case(pass)
      case(0)
         f1 = dt*(32.0/60.0)
         f2 = f1
         f3 = 0.0
         use_local = .false.
         copy_back = .false.
      case(1)
         f1 = dt*( 8.0/60.0)
         f2 = dt*(25.0/60.0)
         f3 = dt*(17.0/60.0)
         use_local = .true.
         copy_back = .true.
      case(2)
         f1 = dt*(20.0/60.0)
         f2 = dt*(45.0/60.0)
         f3 = dt*(25.0/60.0)
         use_local = .true.
         copy_back = .false.
      case default
         return
      end select

      ! Pressure-gradient contributions
      if (.not. use_local) then
         call momxp(durk,pnew)
         call momyp(dvrk,pnew)
         call momzp(dwrk,pnew)
      else
         call momxp(dnewu,pnew)
         call momyp(dnewv,pnew)
         call momzp(dneww,pnew)
      end if

      !$omp parallel default(none)    &
      !$omp&shared(f1,use_local)      &
      !$omp&shared(durk,dvrk,dwrk,dcrk) &
      !$omp&shared(dnewu,dnewv,dneww,dnewc) &
      !$omp&shared(unew ,vnew ,wnew ,cnew ) &
      !$omp&shared(dudt ,dvdt ,dwdt ,dcdt) &
      !$omp&private(i,j,k)
      !$omp do
      do k=1,kmax
          do j=1,jmax
               do i=1,imax
                   if (.not. use_local) then
                      dudt(i,j,k) = unew(i,j,k) + f1*durk(i,j,k)
                      dvdt(i,j,k) = vnew(i,j,k) + f1*dvrk(i,j,k)
                      dwdt(i,j,k) = wnew(i,j,k) + f1*dwrk(i,j,k)
                      dcdt(i,j,k) = cnew(i,j,k) ! + f1*dcrk(i,j,k) ! (SS)
                   else
                      dudt(i,j,k) = unew(i,j,k) + f1*dnewu(i,j,k)
                      dvdt(i,j,k) = vnew(i,j,k) + f1*dnewv(i,j,k)
                      dwdt(i,j,k) = wnew(i,j,k) + f1*dneww(i,j,k)
                      dcdt(i,j,k) = cnew(i,j,k) ! + f1*dnewc(i,j,k) ! (SS)
                   end if
               enddo
          enddo
      enddo
      !$omp end parallel

      ! Advection-diffusion + buoyancy for momentum and scalar advection-diffusion
      if (.not. use_local) then
         call momxad(durk,unew,vnew,wnew)
         call momyad(dvrk,unew,vnew,wnew)
         call momzad(dwrk,unew,vnew,wnew)
         dwrk_old(:,:,:) = dwrk(:,:,:) ! (SS)
         call momzbuoy(dwrk,cnew)
         dwrk(:,:,:) = dwrk(:,:,:) + dwrk_old(:,:,:) ! (SS)
         call scalad(dcrk,unew,vnew,wnew,cnew)
         if (runmode == RUNMODE_TEST_SCALAR_SCHEME) call tests_scalar_scheme_src(dcrk,dwrk,time)
      else
         call momxad(dnewu,unew,vnew,wnew)
         call momyad(dnewv,unew,vnew,wnew)
         call momzad(dneww,unew,vnew,wnew)
         dneww_old(:,:,:) = dneww(:,:,:) ! (SS)
         call momzbuoy(dneww,cnew)
         dneww(:,:,:) = dneww(:,:,:) + dneww_old(:,:,:) ! (SS)
         call scalad(dnewc,unew,vnew,wnew,cnew)
         if (runmode == RUNMODE_TEST_SCALAR_SCHEME) call tests_scalar_scheme_src(dnewc,dneww,time)
      end if

      !$omp parallel default(none)      &
      !$omp&shared(f2,f3,use_local,copy_back) &
      !$omp&shared(durk,dvrk,dwrk,dcrk) &
      !$omp&shared(dnewu,dnewv,dneww,dnewc) &
      !$omp&shared(dudt ,dvdt ,dwdt ,dcdt) &
      !$omp&private(i,j,k)
      !$omp do
      do k=1,kmax
          do j=1,jmax
               do i=1,imax
                   if (.not. use_local) then
                      dudt(i,j,k) = dudt(i,j,k) + f2*durk(i,j,k)
                      dvdt(i,j,k) = dvdt(i,j,k) + f2*dvrk(i,j,k)
                      dwdt(i,j,k) = dwdt(i,j,k) + f2*dwrk(i,j,k)
                      dcdt(i,j,k) = dcdt(i,j,k) + f2*dcrk(i,j,k)
                   else
                      dudt(i,j,k) = dudt(i,j,k) + ( f2*dnewu(i,j,k) - f3*durk(i,j,k) )
                      dvdt(i,j,k) = dvdt(i,j,k) + ( f2*dnewv(i,j,k) - f3*dvrk(i,j,k) )
                      dwdt(i,j,k) = dwdt(i,j,k) + ( f2*dneww(i,j,k) - f3*dwrk(i,j,k) )
                      dcdt(i,j,k) = dcdt(i,j,k) + ( f2*dnewc(i,j,k) - f3*dcrk(i,j,k) )
                      if (copy_back) then
                         durk(i,j,k) = dnewu(i,j,k)
                         dvrk(i,j,k) = dnewv(i,j,k)
                         dwrk(i,j,k) = dneww(i,j,k)
                         dcrk(i,j,k) = dnewc(i,j,k)
                      end if
                   end if
               enddo
          enddo
      enddo
      !$omp end parallel

      ! Wall shear bookkeeping (stage-form specific)
      dVeul_loc = (1.0/dxi)*(1.0/dyi)*(1.0/dzi)
      wallshearnew = 0.0
!!$omp parallel default(shared) &
!!$omp&private(i,j,k) &
!!$omp&reduction(+:wallshearnew)
!!$omp do
      do j=1,jmax
         do i=1,imax
            wallshearnew = wallshearnew + visc*(vnew(i,j,kmax+1)-vnew(i,j,kmax))*dzi*(1.0/dxi)*(1.0/dyi) &
                                                      - visc*(vnew(i,j,1)-vnew(i,j,0))*dzi*(1.0/dxi)*(1.0/dyi)
         enddo
      enddo
!!$omp end parallel
      call mpi_allreduce(wallshearnew,wallshear_all,1,mpi_real8,mpi_sum,comm_cart,error)
      avgShear = wallshear_all/(dVeul_loc*itot*jtot*kmax)
      select case(pass)
      case(0)
         wallshearnew = (f1/dt)*avgShear
         wallshearold = avgShear
      case(1)
         wallshearnew = (f2/dt)*avgShear - (f3/dt)*wallshearold
         wallshearold = avgShear
      case(2)
         wallshearnew = (f2/dt)*avgShear - (f3/dt)*wallshearold
      end select

   end subroutine rk_stage
   ! -----------------------------------------------------------------------
end module mod_rk
