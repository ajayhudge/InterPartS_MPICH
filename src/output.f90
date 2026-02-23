module mod_output
  use mod_common
  use mod_common_mpi
  use decomp_2d
  use decomp_2d_io
  use mod_post
  use mod_phase_indicator
  implicit none
  private
  public post1d,post2d,post3d,historypart
contains
  !
  subroutine post1d(time,istep,nstep)
    implicit none
    integer, intent(in) :: istep,nstep
    real, intent(in) :: time
    real, dimension(0:k1) :: um  ,vm  ,wm   ,pm,  &
                             u2  ,v2  ,w2   ,p2,  &
                             uv  ,uw  ,vw   ,     &
                             upm ,vpm ,wpm  ,     &
                             up2 ,vp2 ,wp2  ,     &
                             upvp,upwp,vpwp ,     &
                             gum ,gvm ,gwm  ,gpm, &
                             gu2 ,gv2 ,gw2  ,gp2, &
                             guv ,guw ,gvw, &
                             ss,ssp, &
                             q1,q2,q3,q4,g1,g2,g3,g4, &
                             cm, c2, uc, vc, wc, cdflx
    real, dimension(0:k1) :: sum_all
    !
    integer :: i,j,k
    real :: dpdx, dpdx_all
    real :: wc_avg, cdflx_avg,Nu,Re_tau,u_tau, Re_RBC,grad_k_wall
    real :: vfluctu,wfluctu,gamvw 
    real, dimension(0:i1,0:j1,0:k1) :: gamu,gamv,gamw,gamp, &
                                       upar,vpar,wpar
    character(len=7) :: istepchar
    !
    if (np>0) then
    call phase_indicator(gamu,gamv,gamw,gamp, &
         upar,vpar,wpar,1)
    else
      gamu = 0
      gamv = 0
      gamw = 0
      gamp = 0
      upar = 0
      vpar = 0
      wpar = 0
    endif

   !  write(istepchar,'(i7.7)') istep
   !  open(unit=200, file='gamw_field_t'//trim(istepchar)//'.dat', status='replace', action='write')
   !  do k = 1, kmax
   !    do j = 1, jmax
   !      write(200, '(1000E12.4)') (gamw(i,j,k), i=1, imax)
   !    end do
   !  end do
   !  close(200)

   !  write(istepchar,'(i7.7)') istep
   !  open(unit=199, file='pressure_field_t'//trim(istepchar)//'.dat', status='replace', action='write')
   !  do k = 1, kmax
   !    do j = 1, jmax
   !      write(199, '(1000E12.4)') (pnew(i,j,k), i=1, imax)
   !    end do
   !  end do
   !  close(199)

   !  write(istepchar,'(i7.7)') istep
   !  open(unit=200, file='gamp_field_t'//trim(istepchar)//'.dat', status='replace', action='write')
   !  do k = 1, kmax
   !    do j = 1, jmax
   !      write(200, '(1000E12.4)') (gamp(i,j,k), i=1, imax)
   !    end do
   !  end do
   !  close(200)


   ! end if

    !$omp parallel default(none) &
    !$omp&private(i,j,k)&
    !$omp&shared(unew,vnew,wnew,pnew,gamu,gamv,gamw,gamp) &
    !$omp&shared(upar,vpar,wpar) &
    !$omp&reduction(+:um,vm,wm,pm,upm,vpm,wpm,gum,gvm,gwm,gpm)
    !$omp do
    do k=1,kmax
       um(k) = 0.
       vm(k) = 0.
       wm(k) = 0.
       pm(k) = 0.
       cm(k) = 0.
       gum(k) = 0.
       gvm(k) = 0.
       gwm(k) = 0.
       gpm(k) = 0.
       upm(k) = 0.
       vpm(k) = 0.
       wpm(k) = 0.
       do j=1,jmax
          do i=1,imax
             if (np > 0) then
             um(k)  = um(k)  + unew(i,j,k)*gamu(i,j,k)
             vm(k)  = vm(k)  + vnew(i,j,k)*gamv(i,j,k)
             wm(k)  = wm(k)  + wnew(i,j,k)*gamw(i,j,k)
             pm(k)  = pm(k)  + pnew(i,j,k)*gamp(i,j,k)
             cm(k)  = cm(k)  + cnew(i,j,k)*gamp(i,j,k) ! SS
             else
             um(k)  = um(k)  + unew(i,j,k)*(1.-gamu(i,j,k))
             vm(k)  = vm(k)  + vnew(i,j,k)*(1.-gamv(i,j,k))
             wm(k)  = wm(k)  + wnew(i,j,k)*(1.-gamw(i,j,k))
             pm(k)  = pm(k)  + pnew(i,j,k)*(1.-gamp(i,j,k))
             cm(k)  = cm(k)  + cnew(i,j,k)*(1.-gamp(i,j,k)) ! SS
             end if             
             upm(k) = upm(k) + upar(i,j,k)*(1.-gamu(i,j,k))
             vpm(k) = vpm(k) + vpar(i,j,k)*(1.-gamv(i,j,k))
             wpm(k) = wpm(k) + wpar(i,j,k)*(1.-gamw(i,j,k))
             gum(k) = gum(k) + gamu(i,j,k)
             gvm(k) = gvm(k) + gamv(i,j,k)
             gwm(k) = gwm(k) + gamw(i,j,k)
             gpm(k) = gpm(k) + gamp(i,j,k)  
          enddo
       enddo
    enddo
    !$omp end parallel
    !
    call mpi_allreduce(um(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    um(:) = sum_all(:)
    !$omp end workshare
    call mpi_allreduce(vm(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    vm(:) = sum_all(:)
    !$omp end workshare
    call mpi_allreduce(wm(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    wm(:) = sum_all(:)
    !$omp end workshare
    call mpi_allreduce(pm(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    pm(:) = sum_all(:)
    !$omp end workshare
     call mpi_allreduce(cm(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    cm(:) = sum_all(:)
    !$omp end workshare
    !
    call mpi_allreduce(upm(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    upm(:) = sum_all(:)
    !$omp end workshare
    call mpi_allreduce(vpm(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    vpm(:) = sum_all(:)
    !$omp end workshare
    call mpi_allreduce(wpm(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    wpm(:) = sum_all(:)
    !$omp end workshare
    !
    call mpi_allreduce(gum(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    gum(:) = sum_all(:)
    !$omp end workshare
    call mpi_allreduce(gvm(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    gvm(:) = sum_all(:)
    !$omp end workshare
    call mpi_allreduce(gwm(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    gwm(:) = sum_all(:)
    !$omp end workshare
    call mpi_allreduce(gpm(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    gpm(:) = sum_all(:)
    !$omp end workshare
    !
    um(0)     = -um(1)
    um(k1)    = -um(kmax)
    vm(0)     = -vm(1)
    vm(k1)    = -vm(kmax)
    wm(0)     = 0.
    wm(kmax)  = 0.
    wm(k1)    = wm(kmax-1)
    pm(0)     = pm(1)
    pm(k1)    = pm(kmax)
    upm(0)    = -upm(1)
    upm(k1)   = -upm(kmax)
    vpm(0)    = -vpm(1)
    vpm(k1)   = -vpm(kmax)
    wpm(0)    = 0.
    wpm(kmax) = 0.
    wpm(k1)   = wm(kmax-1)
    gum(0)    = gum(1)
    gum(k1)   = gum(kmax)
    gvm(0)    = gvm(1)
    gvm(k1)   = gvm(kmax)
    gwm(0)    = gwm(1)
    gwm(k1)   = gwm(kmax-1)
    !
    !$omp parallel default(none) &
    !$omp&private(i,j,k) &
    !$omp&shared(unew,vnew,wnew,pnew,gamu,gamv,gamw,gamp) &
    !$omp&shared(upar,vpar,wpar) &
    !$omp&reduction(+:u2,v2,w2,p2,up2,vp2,wp2,gp2)
    !$omp do
    do k=1,kmax
       u2(k)  = 0.
       v2(k)  = 0.
       w2(k)  = 0.
       up2(k) = 0.
       vp2(k) = 0.
       wp2(k) = 0.
       p2(k)  = 0.
       c2(k)  = 0.
       gp2(k) = 0.
       do j=1,jmax
          do i=1,imax
             if (np > 0) then
             u2(k)  = u2(k) + unew(i,j,k)**2.*gamu(i,j,k)
             v2(k)  = v2(k) + vnew(i,j,k)**2.*gamv(i,j,k)
             w2(k)  = w2(k) + wnew(i,j,k)**2.*gamw(i,j,k)
             p2(k)  = p2(k) + pnew(i,j,k)**2.*gamp(i,j,k)
             c2(k)  = c2(k) + cnew(i,j,k)**2.*gamp(i,j,k) ! SS
             else
             u2(k)  = u2(k) + unew(i,j,k)**2.*(1.-gamu(i,j,k))
             v2(k)  = v2(k) + vnew(i,j,k)**2.*(1.-gamv(i,j,k))
             w2(k)  = w2(k) + wnew(i,j,k)**2.*(1.-gamw(i,j,k))
             p2(k)  = p2(k) + pnew(i,j,k)**2.*(1.-gamp(i,j,k))
             c2(k)  = c2(k) + cnew(i,j,k)**2.*(1.-gamp(i,j,k)) ! SS
             end if
             up2(k) = up2(k) + upar(i,j,k)**2.*(1.-gamu(i,j,k))
             vp2(k) = vp2(k) + vpar(i,j,k)**2.*(1.-gamv(i,j,k))
             wp2(k) = wp2(k) + wpar(i,j,k)**2.*(1.-gamw(i,j,k))
             gp2(k) = gp2(k) + gamp(i,j,k)**2.
          enddo
       enddo
    enddo
    !$omp end parallel
    !
    call mpi_allreduce(u2(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    u2(:) = sum_all(:)
    !$omp end workshare
    call mpi_allreduce(v2(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    v2(:) = sum_all(:)
    !$omp end workshare
    call mpi_allreduce(w2(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    w2(:) = sum_all(:)
    !$omp end workshare
    call mpi_allreduce(p2(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    p2(:) = sum_all(:)
    !$omp end workshare
    call mpi_allreduce(c2(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    c2(:) = sum_all(:)
    !$omp end workshare
    !
    call mpi_allreduce(up2(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    up2(:) = sum_all(:)
    !$omp end workshare
    call mpi_allreduce(vp2(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    vp2(:) = sum_all(:)
    !$omp end workshare
    call mpi_allreduce(wp2(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    wp2(:) = sum_all(:)
    !$omp end workshare
    !
    u2(0)     = -u2(1)
    u2(k1)    = -u2(kmax)
    v2(0)     = -v2(1)
    v2(k1)    = -v2(kmax)
    w2(0)     = 0.
    w2(kmax)  = 0.
    w2(k1)    = w2(kmax-1)
    p2(0)     = p2(1)
    p2(k1)    = p2(kmax)
    up2(0)    = -up2(1)
    up2(k1)   = -up2(kmax)
    vp2(0)    = -vp2(1)
    vp2(k1)   = -vp2(kmax)
    wp2(0)    = 0.
    wp2(kmax) = 0.
    wp2(k1)   = w2(kmax-1)
    gp2(0)    = gp2(1)
    gp2(k1)   = gp2(kmax)
    !
    !$omp parallel default(none) &
    !$omp&private(i,j,k) &
    !$omp&shared(unew,vnew,wnew,gamu,gamv,gamw,gamp) &
    !$omp&shared(upar,vpar,wpar) &
    !$omp&reduction(+:uv,uw,vw,upvp,vpwp,upwp,guv,guw,gvw)
    !$omp do
    do k=1,kmax
       uv(k)=0.
       uw(k)=0.
       vw(k)=0.
       upvp(k)=0.
       upwp(k)=0.
       vpwp(k)=0.
       uc(k)=0.
       vc(k)=0.
       wc(k)=0.
       cdflx(k)=0.
       guv(k)=0.
       guw(k)=0.
       gvw(k)=0.
       do j=1,jmax
          do i=1,imax
             uv(k) = uv(k)     + 0.25*(unew(i,j+1,k) + unew(i,j,k))* &
                                      (vnew(i+1,j,k) + vnew(i,j,k))* ( &
                                 0.25*(gamu(i,j+1,k) + gamu(i,j,k) + &
                                       gamv(i+1,j,k) + gamv(i,j,k))  )
             uw(k) = uw(k)     + 0.25*(unew(i,j,k+1) + unew(i,j,k))* &
                                      (wnew(i+1,j,k) + wnew(i,j,k))* ( &
                                 0.25*(gamu(i,j,k+1) + gamu(i,j,k) + &
                                       gamw(i+1,j,k) + gamw(i,j,k))  )
             vw(k) = vw(k)     + 0.25*(vnew(i,j,k+1) + vnew(i,j,k))* &
                                      (wnew(i,j+1,k) + wnew(i,j,k))* ( &
                                 0.25*(gamv(i,j,k+1) + gamv(i,j,k) + &
                                       gamw(i,j+1,k) + gamw(i,j,k))  )
             
             wc(k) = wc(k)     + 0.5* (cnew(i,j,k+1) + cnew(i,j,k))* &
                                      wnew(i,j,k)* (1 - gamw(i,j,k))  ! SS
             cdflx(k) = cdflx(k) - kappa*(cnew(i,j,k+1) - cnew(i,j,k))*dzi

             upvp(k) = upvp(k) + 0.25*(upar(i,j+1,k) + upar(i,j,k))* &
                                      (vpar(i+1,j,k) + vpar(i,j,k))* ( &
                                 0.25*(2.-gamu(i,j+1,k) - gamu(i,j,k) + &
                                       2.-gamv(i+1,j,k) - gamv(i,j,k))  )
             upwp(k) = upwp(k) + 0.25*(upar(i,j,k+1) + upar(i,j,k))* &
                                      (wpar(i+1,j,k) + wpar(i,j,k))* ( &
                                 0.25*(2.-gamu(i,j,k+1) - gamu(i,j,k) + &
                                       2.-gamw(i+1,j,k) - gamw(i,j,k))  )
             vpwp(k) = vpwp(k) + 0.25*(vpar(i,j,k+1) + vpar(i,j,k))* &
                                      (wpar(i,j+1,k) + wpar(i,j,k))* ( &
                                 0.25*(2.-gamv(i,j,k+1) - gamv(i,j,k) + &
                                       2.-gamw(i,j+1,k) - gamw(i,j,k))  )
             guv(k) = guv(k)   + 0.25*(gamu(i,j+1,k) + gamu(i,j,k) + &
                                      (gamv(i+1,j,k) + gamv(i,j,k)))
             guw(k) = guw(k)   + 0.25*(gamu(i,j,k+1) + gamu(i,j,k) + &
                                      (gamw(i+1,j,k) + gamw(i,j,k)))
             gvw(k) = gvw(k)   + 0.25*(gamv(i,j,k+1) + gamv(i,j,k) + &
                                      (gamw(i,j+1,k) + gamw(i,j,k)))
          enddo
       enddo
    enddo
    !$omp end parallel
    !
    call mpi_allreduce(uv(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    uv(:) = sum_all(:)
    !$omp end workshare
    call mpi_allreduce(uw(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    uw(:) = sum_all(:)
    !$omp end workshare
    call mpi_allreduce(vw(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    vw(:) = sum_all(:)
    !$omp end workshare
    call mpi_allreduce(wc(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    wc(:) = sum_all(:)
    !$omp end workshare
    call mpi_allreduce(cdflx(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    cdflx(:) = sum_all(:)
    !$omp end workshare
    !
    call mpi_allreduce(upvp(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    upvp(:) = sum_all(:)
    !$omp end workshare
    call mpi_allreduce(upwp(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    upwp(:) = sum_all(:)
    !$omp end workshare
    call mpi_allreduce(vpwp(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    vpwp(:) = sum_all(:)
    !$omp end workshare
    !
    call mpi_allreduce(guv(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    guv(:) = sum_all(:)
    !$omp end workshare
    call mpi_allreduce(guw(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    guw(:) = sum_all(:)
    !$omp end workshare
    call mpi_allreduce(gvw(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    gvw(:) = sum_all(:)
    !$omp end workshare
    !
    uv(0)      = -uv(1)
    uv(k1)     = -uv(kmax)
    uw(0)      = 0.
    uw(kmax)   = 0.
    uw(k1)     = uw(kmax-1)
    vw(0)      = 0.
    vw(k1)     = vw(kmax-1)
    vw(kmax)   = 0.
    upvp(0)    = -upvp(1)
    upvp(k1)   = -upvp(kmax)
    upwp(0)    = 0.
    upwp(kmax) = 0.
    upwp(k1)   = upwp(kmax-1)
    vpwp(0)    = 0.
    vpwp(k1)   = vpwp(kmax-1)
    vpwp(kmax) = 0.
    guv(0)     = guv(1)
    guv(k1)    = guv(kmax)
    guw(0)     = guw(1)
    guw(k1)    = guw(kmax)
    gvw(0)     = gvw(1)
    gvw(k1)    = gvw(kmax)
    !
    ! compute mean and rms velocities, Reynolds and viscous stresses
    ! for the fluid and particles, and output them
    if (np > 0) then
    um(1:kmax) = um(1:kmax)/gum(1:kmax)
    vm(1:kmax) = vm(1:kmax)/gvm(1:kmax)
    wm(1:kmax) = wm(1:kmax)/gwm(1:kmax)
    cm(1:kmax) = cm(1:kmax)/gpm(1:kmax) ! SS
    pm(1:kmax) = pm(1:kmax)/gpm(1:kmax)
    u2(1:kmax) = sqrt(u2(1:kmax)/gum(1:kmax)-um(1:kmax)**2.)
    v2(1:kmax) = sqrt(v2(1:kmax)/gvm(1:kmax)-vm(1:kmax)**2.)
    w2(1:kmax) = sqrt(w2(1:kmax)/gwm(1:kmax)-wm(1:kmax)**2.)
    p2(1:kmax) = sqrt(p2(1:kmax)/gpm(1:kmax)-pm(1:kmax)**2.)
    c2(1:kmax) = sqrt(c2(1:kmax)/gpm(1:kmax)-cm(1:kmax)**2.)
    vw(1:kmax) = vw(1:kmax)/gvw(1:kmax)-0.5*(vm(1:kmax)+vm(1+1:kmax+1))*wm(1:kmax)
    vw(0)  = 0.
    else
    um(1:kmax) = um(1:kmax)/(1.*itot*jtot-gum(1:kmax))
    vm(1:kmax) = vm(1:kmax)/(1.*itot*jtot-gvm(1:kmax))
    wm(1:kmax) = wm(1:kmax)/(1.*itot*jtot-gwm(1:kmax))
    cm(1:kmax) = cm(1:kmax)/(1.*itot*jtot-gpm(1:kmax)) ! SS
    pm(1:kmax) = pm(1:kmax)/(1.*itot*jtot-gpm(1:kmax))
    u2(1:kmax) = sqrt(u2(1:kmax)/(1.*itot*jtot-gum(1:kmax))-um(1:kmax)**2.)
    v2(1:kmax) = sqrt(v2(1:kmax)/(1.*itot*jtot-gvm(1:kmax))-vm(1:kmax)**2.)
    w2(1:kmax) = sqrt(w2(1:kmax)/(1.*itot*jtot-gwm(1:kmax))-wm(1:kmax)**2.)
    p2(1:kmax) = sqrt(p2(1:kmax)/(1.*itot*jtot-gpm(1:kmax))-pm(1:kmax)**2.)
    c2(1:kmax) = sqrt(c2(1:kmax)/(1.*itot*jtot-gpm(1:kmax))-cm(1:kmax)**2.)
    vw(1:kmax) = vw(1:kmax)/(1.*itot*jtot-gvw(1:kmax))-0.5*(vm(1:kmax)+vm(1+1:kmax+1))*wm(1:kmax)
    vw(0)  = 0.
    end if

    Re_RBC = maxval(u2)*lz/visc ! RBC Reynolds number (SS)
    wc(1:kmax) = wc(1:kmax)/(1.*itot*jtot-gwm(1:kmax))
    cdflx(1:kmax) = cdflx(1:kmax)/(1.*itot*jtot-gwm(1:kmax))

    vm(0)  = -vm(1)
    vm(k1) = -vm(kmax)
    ss(1:kmax) = (vm(1+1:kmax+1)-vm(1:kmax))*dzi*visc ! Question: why only y-directional shear stress? (SS)
    ss(0) = (vm(1)-vm(0))*dzi*visc
    um(0)  = -um(1)
    um(k1) = -um(kmax)
    ss(1:kmax) = sqrt(ss(1:kmax)**2 + ((um(1+1:kmax+1)-um(1:kmax))*dzi*visc)**2) ! add x-directional shear stress (SS)
    ss(0) = sqrt(ss(0)**2 + ((um(1)-um(0))*dzi*visc)**2)
    grad_k_wall = dzi*(sqrt(0.5*(u2(1)**2 + v2(1)**2 + w2(1)**2)) - sqrt(0.5*(u2(0)**2 + v2(0)**2 + w2(0)**2)))
    ! turbulent kinetic energy root gradient at the wall (SS)
    u_tau = sqrt(visc*(grad_k_wall)) ! shear velocity at the wall (SS)
    Re_tau = u_tau*lz/visc ! friction Reynolds number (SS)   

    ! Calculate the Nusselt number
    wc_avg = sum(wc(1:kmax))/kmax
    cdflx_avg = sum(cdflx(1:kmax))/kmax
    Nu = lz/(kappa*(bc_c_bot_val-bc_c_top_val)) * (wc_avg + cdflx_avg)


    write(istepchar,'(i7.7)') istep
    if (myid .eq. 0) then
       open(unit=30,file='info_1dstat.out')
       write(30,'(8E15.7)') 1.*istep,time,visc,1.*ktot,dx,uref,lref,tref
       close(30)
      !  open(30,file='stats_fluid_uvw_'//istepchar//'.out')
      !  do k=0,k1
      !     write(30,'(12E15.7)') 1.*k, & 
      !          um(k),vm(k),wm(k),pm(k), &
      !          u2(k),v2(k),w2(k),p2(k), & 
      !          uv(k),uw(k),vw(k)
      !  enddo
      !  close(30)
      !  open(30,file='stats_fluid_c_'//istepchar//'.out')
      !  do k=0,k1
      !     write(30,'(12E15.7)') 1.*k, zc(k), & 
      !          cm(k),c2(k),uc(k),vc(k), &
      !          wc(k),cdflx(k)
      !  enddo
      !  close(30)

      ! (SS) 2025/11/11: commented out plotting to avoid excessive output files 
      ! if (np > 0) then       
      !  open(30,file='stats_phase_uvw_'//istepchar//'.out')
      !  do k=0,k1
      !     write(30,'(8E15.7)') 1.*k, & 
      !          gum(k),gvm(k),gwm(k),gpm(k), &
      !          guv(k),guw(k),gvw(k)
      !  enddo
      !  close(30)
      ! end if
      
      ! if (np > 0) then
      !  open(30,file='stats_parts_uvw_'//istepchar//'.out')
      !  do k=0,k1
      !     write(30,'(10E15.7)') 1.*k, & 
      !          upm(k) ,vpm(k) ,wpm(k), &
      !          up2(k) ,vp2(k) ,wp2(k), &
      !          upvp(k),upwp(k),vpwp(k)
      !  enddo
      !  close(30)
      ! end if
      !
       
       !
      !  open(30,file='stats_fluid_uvw_now.out')
      !  do k=1,kmax
      !     write(30,'(13E15.7)') (k-0.5)/dzi, (k-0.5)/dzi/visc, &
      !          um(k) ,vm(k) ,0.5*(wm(k)+wm(k-1)), &
      !          u2(k) ,v2(k) ,0.5*(w2(k)+w2(k-1)), &
      !          pm(k) ,p2(k) , &
      !          0.5*(vw(k)+vw(k-1)) ,0.5*(ss(k)+ss(k-1)), &
      !          0.5*(vw(k)+vw(k-1)) + 0.5*(ss(k)+ss(k-1))
      !  enddo
      !  close(30)

      ! (SS) 2025/11/11: commented out plotting to avoid excessive output files
      !  open(30,file='Nu_vs_time.out',status='unknown',position='append')
      !     write(30,'(10E15.7)') 1.*time, Nu
      !  close(30)

      !  open(30,file='Re_vs_time.out',status='unknown',position='append')
      !     write(30,'(10E15.7)') 1.*time, Re_RBC
      !  close(30)     

      !  open(30,file='Re_tau_vs_time.out',status='unknown',position='append')
      !     write(30,'(10E15.7)') 1.*time, Re_tau
      !  close(30) 
      !

       !$omp workshare
       upm(1:kmax) = upm(1:kmax)/(1.*itot*jtot-gum(1:kmax))
       vpm(1:kmax) = vpm(1:kmax)/(1.*itot*jtot-gvm(1:kmax))
       wpm(1:kmax) = wpm(1:kmax)/(1.*itot*jtot-gwm(1:kmax))
       up2(1:kmax) = sqrt(up2(1:kmax)/(1.*itot*jtot-gum(1:kmax))-upm(1:kmax)**2.)
       vp2(1:kmax) = sqrt(vp2(1:kmax)/(1.*itot*jtot-gvm(1:kmax))-vpm(1:kmax)**2.)
       wp2(1:kmax) = sqrt(wp2(1:kmax)/(1.*itot*jtot-gwm(1:kmax))-wpm(1:kmax)**2.)
       vpwp(1:kmax) = vpwp(1:kmax)/(1.*itot*jtot-gvw(1:kmax)) - &
                      0.5*(vpm(1:kmax)+vpm(1+1:kmax+1))*wpm(1:kmax)
       !$omp end workshare
       vpwp(0) = 0.
       vpm(0)  = -vpm(1)
       vpm(k1) = -vpm(kmax)
       ssp(1:kmax) = (vpm(1+1:kmax+1)-vpm(1:kmax))*dzi*visc
       ssp(0) = (vpm(1)-vpm(0))*dzi*visc
       !
       ! compute mean and rms velocities, Reynolds and viscous stresses
       ! for the fluid and particles, and plot them
       !
      ! (SS) 2025/11/11: commented out plotting to avoid excessive output files
      ! if (np > 0) then
      !  open(30,file='stats_parts_uvw_now.out')
      !  do k=1,kmax
      !     write(30,'(11E15.7)') (k-0.5)/dzi   ,(k-0.5)/dzi/visc, &
      !          upm(k)  ,vpm(k) ,0.5*(wpm(k)+wpm(k-1)), &
      !          up2(k)  ,vp2(k) ,0.5*(wp2(k)+wp2(k-1)), &
      !          0.5*(vpwp(k)+vpwp(k-1)), 0.5*(ssp(k)+ssp(k-1)), &
      !          0.5*(vpwp(k)+vpwp(k-1))+ 0.5*(ssp(k)+ssp(k-1))
      !  enddo
      !  close(30)
      !  open(30,file='stats_phase_uvw_now.out')
      !  do k=1,kmax
      !     write(30,'(3E15.7)') (k-0.5)/dzi, (k-0.5)/dzi/visc,(1.-gpm(k)/(1.*itot*jtot))
      !  enddo
      !  close(30)
      ! end if
      ! 
    endif
    !
    call outputpart(istep) ! partvisu*.bin: output acce, vel, pos, etc of particles (SS)
    !    

    ! no quadrant analysis please
    return
    !
    ! quadrant analysis
    !
    !$omp parallel default(none) &
    !$omp&private(i,j,k,vfluctu,wfluctu,gamvw) &
    !$omp&shared(vnew,wnew,gamv,gamw,vm,wm) &
    !$omp&reduction(+:q1,q2,q3,q4,g1,g2,g3,g4)
    !$omp do
    do k=1,kmax
       q1(k) = 0.
       q2(k) = 0.
       q3(k) = 0.
       q4(k) = 0.
       g1(k) = 0.
       g2(k) = 0.
       g3(k) = 0.
       g4(k) = 0.
       do j=1,jmax
          do i=1,imax
             vfluctu = vnew(i,j,k) - vm(k)
             wfluctu = wnew(i,j,k) - wm(k)
             gamvw   =  0.25*( gamv(i,j,k+1) + gamv(i,j,k) + &
                               gamw(i,j+1,k) + gamw(i,j,k) )
             if    (vfluctu.gt.0.and.wfluctu.gt.0) then ! 1st quadrant
               q1(k) = q1(k) + vfluctu*wfluctu*gamvw
               g1(k) = g1(k) + gamvw
             elseif(vfluctu.lt.0.and.wfluctu.gt.0) then ! 2nd quadrant
               q2(k) = q2(k) + vfluctu*wfluctu*gamvw
               g1(k) = g1(k) + gamvw
             elseif(vfluctu.lt.0.and.wfluctu.lt.0) then ! 3rd quadrant
               q3(k) = q3(k) + vfluctu*wfluctu*gamvw
               g1(k) = g1(k) + gamvw
             elseif(vfluctu.gt.0.and.wfluctu.lt.0) then ! 4th quadrant
               q4(k) = q4(k) + vfluctu*wfluctu*gamvw
               g1(k) = g1(k) + gamvw
             endif
          enddo
       enddo
    enddo
    !$omp end parallel
    call mpi_allreduce(q1(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    q1(:) = sum_all(:)
    !$omp end workshare
    call mpi_allreduce(q2(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    q2(:) = sum_all(:)
    !$omp end workshare
    call mpi_allreduce(q3(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    q3(:) = sum_all(:)
    !$omp end workshare
    call mpi_allreduce(q4(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    q4(:) = sum_all(:)
    !$omp end workshare
    call mpi_allreduce(g1(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    g1(:) = sum_all(:)
    !$omp end workshare
    call mpi_allreduce(g2(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    g2(:) = sum_all(:)
    !$omp end workshare
    call mpi_allreduce(g3(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    g3(:) = sum_all(:)
    !$omp end workshare
    call mpi_allreduce(g4(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    g4(:) = sum_all(:)
    !$omp end workshare
    !
    call phase_indicator(gamu,gamv,gamw,gamp, &
         upar,vpar,wpar,0)
    call vorticity(unew, &
                   vnew, &
                   wnew, &
                   gamu(1:imax,1:jmax,1:kmax), &
                   gamv(1:imax,1:jmax,1:kmax), &
                   gamw(1:imax,1:jmax,1:kmax)) ! gamu,gamv,gamw -> vorticity at cell center
    !$omp parallel default(none) &
    !$omp&private(i,j,k)&
    !$omp&shared(unew,vnew,wnew,pnew,gamu,gamv,gamw,gamp) &
    !$omp&shared(upar,vpar,wpar) &
    !$omp&reduction(+:um,vm,wm,upm,vpm,wpm,gpm)
    !$omp do
    do k=1,kmax
       um(k) = 0.
       vm(k) = 0.
       wm(k) = 0.
       gpm(k) = 0.
       upm(k) = 0.
       vpm(k) = 0.
       wpm(k) = 0.
       do j=1,jmax
          do i=1,imax
             um(k)  = um(k)  + gamu(i,j,k)*gamp(i,j,k)
             vm(k)  = vm(k)  + gamv(i,j,k)*gamp(i,j,k)
             wm(k)  = wm(k)  + gamw(i,j,k)*gamp(i,j,k)
             upm(k) = upm(k) + upar(i,j,k)*(1.-gamp(i,j,k))
             vpm(k) = vpm(k) + vpar(i,j,k)*(1.-gamp(i,j,k))
             wpm(k) = wpm(k) + wpar(i,j,k)*(1.-gamp(i,j,k))
             gpm(k) = gpm(k) + gamp(i,j,k)
          enddo
       enddo
    enddo
    !$omp end parallel
!
    call mpi_allreduce(um(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    um(:) = sum_all(:)
    !$omp end workshare
    call mpi_allreduce(vm(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    vm(:) = sum_all(:)
    !$omp end workshare
    call mpi_allreduce(wm(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    wm(:) = sum_all(:)
    !$omp end workshare
    !
    call mpi_allreduce(upm(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    upm(:) = sum_all(:)
    !$omp end workshare
    call mpi_allreduce(vpm(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    vpm(:) = sum_all(:)
    !$omp end workshare
    call mpi_allreduce(wpm(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    wpm(:) = sum_all(:)
    !$omp end workshare
    !
    call mpi_allreduce(gpm(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    gpm(:) = sum_all(:)
    !$omp end workshare
    !
    !$omp parallel default(none) &
    !$omp&private(i,j,k)&
    !$omp&shared(gamu,gamv,gamw,gamp) &
    !$omp&shared(upar,vpar,wpar) &
    !$omp&reduction(+:u2,v2,w2,up2,vp2,wp2)
    !$omp do
    do k=1,kmax
       u2(k) = 0.
       v2(k) = 0.
       w2(k) = 0.
       up2(k) = 0.
       vp2(k) = 0.
       wp2(k) = 0.
       do j=1,jmax
          do i=1,imax
             u2(k) = u2(k) + gamu(i,j,k)**2.*gamp(i,j,k)
             v2(k) = v2(k) + gamv(i,j,k)**2.*gamp(i,j,k)
             w2(k) = w2(k) + gamw(i,j,k)**2.*gamp(i,j,k)
             up2(k) = up2(k) + upar(i,j,k)**2.*(1.-gamp(i,j,k))
             vp2(k) = vp2(k) + vpar(i,j,k)**2.*(1.-gamp(i,j,k))
             wp2(k) = wp2(k) + wpar(i,j,k)**2.*(1.-gamp(i,j,k))
          enddo
       enddo
    enddo
    !$omp end parallel
    !
    call mpi_allreduce(u2(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    u2(:) = sum_all(:)
    !$omp end workshare
    call mpi_allreduce(v2(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    v2(:) = sum_all(:)
    !$omp end workshare
    call mpi_allreduce(w2(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    w2(:) = sum_all(:)
    !$omp end workshare
    !
    call mpi_allreduce(up2(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    up2(:) = sum_all(:)
    !$omp end workshare
    call mpi_allreduce(vp2(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    vp2(:) = sum_all(:)
    !$omp end workshare
    call mpi_allreduce(wp2(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    wp2(:) = sum_all(:)
    !$omp end workshare
    !
    ! Find lazyness in the loop below 
    !
    !$omp parallel default(none) &
    !$omp&private(i,j,k)&
    !$omp&shared(unew,vnew,wnew,pnew,gamu,gamv,gamw,gamp) &
    !$omp&shared(upar,vpar,wpar) &
    !$omp&reduction(+:uv,uw,vw,upvp,upwp,vpwp)
    !$omp do
    do k=1,kmax
       uv(k)=0.
       uw(k)=0.
       vw(k)=0.
       upvp(k)=0.
       upwp(k)=0.
       vpwp(k)=0.
       do j=1,jmax
          do i=1,imax
             uv(k) = uv(k)     + 0.25*(unew(i,j+0,k) + unew(i,j,k))* &
                                      (vnew(i+0,j,k) + vnew(i,j,k))* ( &
                                 0.25*(gamp(i,j+0,k) + gamp(i,j,k) + &
                                       gamp(i+0,j,k) + gamp(i,j,k))  )
             uw(k) = uw(k)     + 0.25*(unew(i,j,k+0) + unew(i,j,k))* &
                                      (wnew(i+0,j,k) + wnew(i,j,k))* ( &
                                 0.25*(gamp(i,j,k+0) + gamp(i,j,k) + &
                                       gamp(i+0,j,k) + gamp(i,j,k))  )
             vw(k) = vw(k)     + 0.25*(vnew(i,j,k+0) + vnew(i,j,k))* &
                                      (wnew(i,j+0,k) + wnew(i,j,k))* ( &
                                 0.25*(gamp(i,j,k+0) + gamp(i,j,k) + &
                                       gamp(i,j+0,k) + gamp(i,j,k))  )
             upvp(k) = upvp(k) + 0.25*(upar(i,j+0,k) + upar(i,j,k))* &
                                      (vpar(i+0,j,k) + vpar(i,j,k))* ( &
                                 0.25*(2.-gamp(i,j+0,k) - gamp(i,j,k) + &
                                       2.-gamp(i+0,j,k) - gamp(i,j,k))  )
             upwp(k) = upwp(k) + 0.25*(upar(i,j,k+0) + upar(i,j,k))* &
                                      (wpar(i+0,j,k) + wpar(i,j,k))* ( &
                                 0.25*(2.-gamp(i,j,k+0) - gamp(i,j,k) + &
                                       2.-gamp(i+0,j,k) - gamp(i,j,k))  )
             vpwp(k) = vpwp(k) + 0.25*(vpar(i,j,k+0) + vpar(i,j,k))* &
                                      (wpar(i,j+0,k) + wpar(i,j,k))* ( &
                                 0.25*(2.-gamp(i,j,k+0) - gamp(i,j,k) + &
                                       2.-gamp(i,j+0,k) - gamp(i,j,k))  )
          enddo
       enddo
    enddo
    !$omp end parallel
    !
    call mpi_allreduce(uv(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    uv(:) = sum_all(:)
    !$omp end workshare
    call mpi_allreduce(uw(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    uw(:) = sum_all(:)
    !$omp end workshare
    call mpi_allreduce(vw(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    vw(:) = sum_all(:)
    !$omp end workshare
    !
    call mpi_allreduce(upvp(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    upvp(:) = sum_all(:)
    !$omp end workshare
    call mpi_allreduce(upwp(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    upwp(:) = sum_all(:)
    !$omp end workshare
    call mpi_allreduce(vpwp(0),sum_all(0),k1+1,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    vpwp(:) = sum_all(:)
    !$omp end workshare
    !
    write(istepchar,'(i7.7)') istep
    if (myid .eq. 0) then
       open(30,file='stats_fluid_rot_'//istepchar//'.out')
       do k=1,kmax
          write(30,'(10E15.7)') 1.*k, & 
               um(k),vm(k),wm(k), &
               u2(k),v2(k),w2(k), &
               uv(k),uw(k),vw(k)
       enddo
       close(30)
       !
       open(30,file='stats_part_rot_'//istepchar//'.out')
       do k=1,kmax
          write(30,'(10E15.7)') 1.*k, & 
               upm(k) ,vpm(k) ,wpm(k), &
               up2(k) ,vp2(k) ,wp2(k), &
               upvp(k),upwp(k),vpwp(k)
       enddo
       close(30)
       !
       !$omp workshare
       um(1:kmax) = um(1:kmax)/gpm(1:kmax)
       vm(1:kmax) = vm(1:kmax)/gpm(1:kmax)
       wm(1:kmax) = wm(1:kmax)/gpm(1:kmax)
       u2(1:kmax) = sqrt(u2(1:kmax)/gpm(1:kmax)-um(1:kmax)**2.)
       v2(1:kmax) = sqrt(v2(1:kmax)/gpm(1:kmax)-vm(1:kmax)**2.)
       w2(1:kmax) = sqrt(w2(1:kmax)/gpm(1:kmax)-wm(1:kmax)**2.)
       vw(1:kmax) = vw(1:kmax)/gvw(1:kmax)-vm(1:kmax)*wm(1:kmax)
       !$omp end workshare
       open(30,file='stats_fluid_rot_now.out')
       do k=1,kmax
          write(30,'(9E15.7)') (k-0.5)/dzi, (k-0.5)/dzi/visc, &
               um(k) ,vm(k) ,wm(k), &
               u2(k) ,v2(k) ,w2(k), &
               vw(k)
       enddo
       close(30)
       !
       !$omp workshare
       upm(1:kmax) = upm(1:kmax)/(1.*itot*jtot-gpm(1:kmax))
       vpm(1:kmax) = vpm(1:kmax)/(1.*itot*jtot-gpm(1:kmax))
       wpm(1:kmax) = wpm(1:kmax)/(1.*itot*jtot-gpm(1:kmax))
       up2(1:kmax) = sqrt(up2(1:kmax)/(1.*itot*jtot-gpm(1:kmax))-upm(1:kmax)**2.)
       vp2(1:kmax) = sqrt(vp2(1:kmax)/(1.*itot*jtot-gpm(1:kmax))-vpm(1:kmax)**2.)
       wp2(1:kmax) = sqrt(wp2(1:kmax)/(1.*itot*jtot-gpm(1:kmax))-wpm(1:kmax)**2.)
       vpwp(1:kmax) = vpwp(1:kmax)/(1.*itot*jtot-gpm(1:kmax)) - &
                      vpm(1:kmax)*wpm(1:kmax)
       !$omp end workshare
       open(30,file='stats_part_rot_now.out')
       do k=1,kmax
          write(30,'(9E15.7)') (k-0.5)/dzi   ,(k-0.5)/dzi/visc, &
               upm(k)  ,vpm(k) ,wpm(k), &
               up2(k)  ,vp2(k) ,wp2(k), &
               vpwp(k)
       enddo
       close(30)
    endif

    return
  end subroutine post1d
  !
  subroutine post2d(istep)
    implicit none
    integer, intent(in) :: istep
    integer :: i,j,k
    real, dimension(0:i1,0:j1,0:k1) :: gamu,gamv,gamw,gamp, &
         upart,vpart,wpart
    !
    call phase_indicator(gamu,gamv,gamw,gamp, &
         upart,vpart,wpart,1)
    !
    !$omp workshare
    gamu(1:imax,1:jmax,1:kmax) = unew(1:imax,1:jmax,1:kmax)*gamu(1:imax,1:jmax,1:kmax)+upart(1:imax,1:jmax,1:kmax)*(1.-gamu(1:imax,1:jmax,1:kmax))
    gamv(1:imax,1:jmax,1:kmax) = vnew(1:imax,1:jmax,1:kmax)*gamv(1:imax,1:jmax,1:kmax)+vpart(1:imax,1:jmax,1:kmax)*(1.-gamv(1:imax,1:jmax,1:kmax))
    gamw(1:imax,1:jmax,1:kmax) = wnew(1:imax,1:jmax,1:kmax)*gamw(1:imax,1:jmax,1:kmax)+wpart(1:imax,1:jmax,1:kmax)*(1.-gamw(1:imax,1:jmax,1:kmax))
    !$omp end workshare
    call write2dplane(gamu(1:imax,1:jmax,1:kmax),1,itot,'vex',istep) ! 'end' of the channel in the spanwise direction
    call write2dplane(gamv(1:imax,1:jmax,1:kmax),1,itot,'vey',istep) ! 'end' of the channel in the spanwise direction
    call write2dplane(gamw(1:imax,1:jmax,1:kmax),1,itot,'vez',istep) ! 'end' of the channel in the spanwise direction
    call write2dplane(gamu(1:imax,1:jmax,1:kmax),2,1,'vex',istep) ! 'begining' of the channel in the streamwise direction
    call write2dplane(gamv(1:imax,1:jmax,1:kmax),2,1,'vey',istep) ! 'begining' of the channel in the streamwise direction
    call write2dplane(gamw(1:imax,1:jmax,1:kmax),2,1,'vez',istep) ! 'begining' of the channel in the streamwise direction
    call write2dplane(gamu(1:imax,1:jmax,1:kmax),3,nint(radius*2*1.5*dzi),'vex',istep) ! 1.5 diameters away from the wall
    call write2dplane(gamv(1:imax,1:jmax,1:kmax),3,nint(radius*2*1.5*dzi),'vey',istep) ! 1.5 diameters away from the wall
    call write2dplane(gamw(1:imax,1:jmax,1:kmax),3,nint(radius*2*1.5*dzi),'vez',istep) ! 1.5 diameters away from the wall
    gamp(1:imax,1:jmax,1:kmax) = gamp(1:imax,1:jmax,1:kmax)*(1.-gamp(1:imax,1:jmax,1:kmax))
    call write2dplane(gamp(1:imax,1:jmax,1:kmax),1,itot,'gmp',istep) ! 'end' of the channel in the spanwise direction
    call write2dplane(gamp(1:imax,1:jmax,1:kmax),2,1,'gmp',istep) ! 'begining' of the channel in the streamwise direction
    call write2dplane(gamp(1:imax,1:jmax,1:kmax),3,nint(radius*2*1.5*dzi),'gmp',istep) ! 1.5 diameters away from the wall
    !
    call phase_indicator(gamu,gamv,gamw,gamp, &
         upart,vpart,wpart,0)
    call vorticity(unew, &
         vnew, &
         wnew, &
         gamu(1:imax,1:jmax,1:kmax), &
         gamv(1:imax,1:jmax,1:kmax), &
         gamw(1:imax,1:jmax,1:kmax)) ! gamu,gamv,gamw -> vorticity at cell center
    !
    !$omp workshare
    gamu(1:imax,1:jmax,1:kmax) = gamu(1:imax,1:jmax,1:kmax)*gamp(1:imax,1:jmax,1:kmax)+upart(1:imax,1:jmax,1:kmax)*(1.-gamp(1:imax,1:jmax,1:kmax))
    gamv(1:imax,1:jmax,1:kmax) = gamv(1:imax,1:jmax,1:kmax)*gamp(1:imax,1:jmax,1:kmax)+vpart(1:imax,1:jmax,1:kmax)*(1.-gamp(1:imax,1:jmax,1:kmax))
    gamw(1:imax,1:jmax,1:kmax) = gamw(1:imax,1:jmax,1:kmax)*gamp(1:imax,1:jmax,1:kmax)+wpart(1:imax,1:jmax,1:kmax)*(1.-gamp(1:imax,1:jmax,1:kmax))
    !$omp end workshare
    call write2dplane(gamu(1:imax,1:jmax,1:kmax),1,itot,'vox',istep) ! 'end' of the channel in the spanwise direction
    call write2dplane(gamv(1:imax,1:jmax,1:kmax),1,itot,'voy',istep) ! 'end' of the channel in the spanwise direction
    call write2dplane(gamw(1:imax,1:jmax,1:kmax),1,itot,'voz',istep) ! 'end' of the channel in the spanwise direction
    call write2dplane(gamu(1:imax,1:jmax,1:kmax),2,1,'vox',istep) ! 'begining' of the channel in the streamwise direction
    call write2dplane(gamv(1:imax,1:jmax,1:kmax),2,1,'voy',istep) ! 'begining' of the channel in the streamwise direction
    call write2dplane(gamw(1:imax,1:jmax,1:kmax),2,1,'voz',istep) ! 'begining' of the channel in the streamwise direction
    call write2dplane(gamu(1:imax,1:jmax,1:kmax),3,nint(radius*2*1.5*dzi),'vox',istep) ! 1.5 diameters away from the wall
    call write2dplane(gamv(1:imax,1:jmax,1:kmax),3,nint(radius*2*1.5*dzi),'voy',istep) ! 1.5 diameters away from the wall
    call write2dplane(gamw(1:imax,1:jmax,1:kmax),3,nint(radius*2*1.5*dzi),'voz',istep) ! 1.5 diameters away from the wall
    !
    return
  end subroutine post2d
  !
  subroutine write2dplane(var,norm,islice,name,istep)
    implicit none
    real, intent(in), dimension(imax,jmax,kmax) :: var
    integer, intent(in) :: norm,islice,istep
    character(len=3) :: name
    character(len=4) :: slicechar
    character(len=7) :: fldnum

    write(fldnum,'(i7.7)') istep
    write(slicechar,'(i4.4)') islice
    select case(norm)
    case(1) !normal to x --> yz plane
       call decomp_2d_write_plane(3,var,norm,islice,name//'_xeq_'//slicechar//'_fld_'//fldnum//'_2d.bin')
    case(2) !normal to y --> zx plane
       call decomp_2d_write_plane(3,var,norm,islice,name//'_yeq_'//slicechar//'_fld_'//fldnum//'_2d.bin')
    case(3) !normal to z --> xy plane
       call decomp_2d_write_plane(3,var,norm,islice,name//'_zeq_'//slicechar//'_fld_'//fldnum//'_2d.bin')
    end select

    return
  end subroutine write2dplane
  !
  subroutine post3d(istep)
    implicit none
    integer, intent(in) :: istep
    real, allocatable, dimension(:,:,:) :: var1,var2,var3,var4
    integer, parameter :: nprocs=dims(1)*dims(2),ksol=kmax/nprocs
    integer :: p
   !  integer :: i,j,k
   !  character(len=7) :: istepchar
   !  real, dimension(0:i1,0:j1,0:k1) :: wpar    
    allocate(var1(1:imax,1:jmax,1:kmax))
    allocate(var3(1:imax,1:jmax,1:kmax))
    allocate(var4(1:imax,1:jmax,1:kmax))
    !
    ! pressure
    !    
    !$omp workshare
    var1(:,:,:) = pnew(1:imax,1:jmax,1:kmax)
    !$omp end workshare
    call write3dscal(istep,imax,jmax,kmax,var1,'pre',1,1,1)
   !  !
   !  ! u velocity
   !  !
   !  !$omp workshare
   !  var1(:,:,:) = 0.5*(unew(1:imax,1:jmax,1:kmax)+unew(1-1:imax-1,1:jmax,1:kmax))
   !  !$omp end workshare
   !  call write3dscal(istep,imax,jmax,kmax,var1,'vex',1,1,1)
   !  !
   !  ! v velocity
   !  !
   !  !$omp workshare
   !  var1(:,:,:) = 0.5*(vnew(1:imax,1:jmax,1:kmax)+vnew(1:imax,1-1:jmax-1,1:kmax))
   !  !$omp end workshare
   !  call write3dscal(istep,imax,jmax,kmax,var1,'vey',1,1,1)
   !  !
    !
    ! w velocity
    !
    !$omp workshare
    var1(:,:,:) = 0.5*(wnew(1:imax,1:jmax,1:kmax)+wnew(1:imax,1:jmax,1-1:kmax-1))
    !$omp end workshare
    call write3dscal(istep,imax,jmax,kmax,var1,'vez',1,1,1) ! 4,4,4
    !
    ! scalar concentration
    !
    !$omp workshare
    var1(:,:,:) = 0.5*(cnew(1:imax,1:jmax,1:kmax)+cnew(1:imax,1:jmax,1-1:kmax-1))
    !$omp end workshare
    call write3dscal(istep,imax,jmax,kmax,var1,'con',1,1,1)
    !
    !
    ! vorticity (SS try)
    !
    call vorticity(unew,vnew,wnew,var1,var3,var4)
    call write3dscal(istep,imax,jmax,kmax,var1,'vox',1,1,1)
    call write3dscal(istep,imax,jmax,kmax,var3,'voy',1,1,1)
   !  call write3dscal(istep,imax,jmax,kmax,var4,'voz',1,1,1)
   !  !
   !  ! dissipation (rate-of-strain)
   !  !
   !  call strain_rate(unew,vnew,wnew,var1)
   !  call write3dscal(istep,imax,jmax,kmax,var1,'str',1,1,1)
   !  !
   !  ! enstrophy
   !  !
   !  call enstrophy(unew,vnew,wnew,var1)
   !  call write3dscal(istep,imax,jmax,kmax,var1,'ens',1,1,1)
   !  !
    deallocate(var1)
    deallocate(var3)
    deallocate(var4)
    !
    allocate(var1(0:i1,0:j1,0:k1),var2(0:i1,0:j1,0:k1))
    call phase_indicator(var2,var2,var2,var1, &
         var2,var2,var2,1)
    !
    ! particles' phase indicator
    !
    call write3dscal(istep,imax,jmax,kmax,var1(1:imax,1:jmax,1:kmax),'sol',1,1,1)
    !call write3dscal(1,imax,jmax,kmax,unew(1:imax,1:jmax,1:kmax),'usm',4,4,4)
    !call write3dscal(1,imax,jmax,kmax,vnew(1:imax,1:jmax,1:kmax),'vsm',4,4,4)
    !call write3dscal(1,imax,jmax,kmax,wnew(1:imax,1:jmax,1:kmax),'wsm',4,4,4)
    deallocate(var1,var2)
    !!
    !! Q-criterion
    !!
    !allocate(var3(1:imax,1:jmax,1:kmax))
    !call q_criterion(var2,var1,var3)
    !call write3dscal(istep,imax,jmax,kmax,var3,'qcr')
    !!
    !! R (third invariant of the velocity-gradient tensor
    !!
    !allocate(var4(1:imax,1:jmax,1:kmax))
    !call compute_r(unew,vnew,wnew,var4)
    !!
    !! swirling strength (lambda_{ci})
    !!
    !call swirl(var3,var4,var1)
    !call write3dscal(istep,imax,jmax,kmax,var3,'swr')
    !!
    return
  end subroutine post3d
  !
  subroutine write3dscal(istep,n1,n2,n3,var,name,iskip,jskip,kskip)
    implicit none
    integer, intent(in) :: istep,n1,n2,n3,iskip,jskip,kskip
    real, intent(in), dimension(n1,n2,n3) :: var
    character(len=3), intent(in) :: name
    integer :: fh
    integer(kind=MPI_OFFSET_KIND) :: filesize,disp
    character :: istepchar*7
    !
    write(istepchar,'(i7.7)') istep
    !call MPI_FILE_OPEN(MPI_COMM_WORLD, datadir//name//istepchar, &
    !     MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, error)
    !filesize = 0_MPI_OFFSET_KIND
    !call MPI_FILE_SET_SIZE(fh,filesize,error)  ! guarantee overwriting
    !disp = 0_MPI_OFFSET_KIND
    !call decomp_2d_write_var(fh,disp,3,var)
    !call MPI_FILE_CLOSE(fh,error)
    call MPI_FILE_OPEN(MPI_COMM_WORLD, name//istepchar//'_3d.bin', &
         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, error)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,error)  ! guarantee overwriting
    disp = 0_MPI_OFFSET_KIND
    call decomp_2d_write_every(3,var,iskip,jskip,kskip,name//istepchar//'_3d.bin',.true.)
    call MPI_FILE_CLOSE(fh,error)
    !
    return
  end subroutine write3dscal
  !
  subroutine outputpart(nr)
    implicit none
    integer :: skip,skipacc
    integer,dimension(0:dims(1)*dims(2)-1) :: npmstr_glob,npmstr_glob_all
    real, allocatable, dimension(:,:) :: posp,velp,ridp,angpos,angvel,accel,angaccel,extravec,extrascal,fltot 
    real :: aux_x,aux_y,aux_z
    integer i,j,p,idp
    integer :: fh
    integer     nr
    real :: xdum
    integer :: lenr
    integer(kind=MPI_OFFSET_KIND) :: filesize,disp
    integer :: mydisp
    character(len=7) :: istepchar
    !
    !inquire (iolength=lenr) xdum
    lenr = sizeof(xdum)
    write(istepchar,'(i7.7)') nr
    !
    ! write particle related data directly in parallel with MPI-IO
    !
    npmstr_glob(:) = 0
    npmstr_glob(myid) = npmstr
    call MPI_ALLREDUCE(npmstr_glob(0),npmstr_glob_all(0),product(dims),MPI_INTEGER,MPI_SUM,comm_cart,error)
    mydisp = 0
    if(myid.ne.0) mydisp = sum(npmstr_glob_all(0:myid-1))
    !allocate(posp(3,npmstr))
    !allocate(velp(3,npmstr))
    !allocate(ridp(1,npmstr))
    !allocate(angpos(3,npmstr))
    !allocate(angvel(3,npmstr))
    allocate(fltot(3,npmstr))
    allocate(accel(3,npmstr))
    !allocate(angaccel(3,npmstr))
    !allocate(extravec(3,npmstr))
    !allocate(extrascal(1,npmstr))
    i = 0
    do p=1,pmax
       if(ap(p)%mslv.gt.0) then
          idp = ap(p)%mslv
          i = i + 1
         !  posp(1,i)      = ap(p)%x
         !  posp(2,i)      = ap(p)%y
         !  posp(3,i)      = ap(p)%z
         !  ridp(1,i)      = 1.*ap(p)%mslv
         !  velp(1,i)      = ap(p)%u 
         !  velp(2,i)      = ap(p)%v
         !  velp(3,i)      = ap(p)%w
         !  aux_x          = sin(ap(p)%theta)*cos(ap(p)%phi)
         !  aux_y          = sin(ap(p)%theta)*sin(ap(p)%phi)
         !  aux_z          = cos(ap(p)%phi)
         !  angpos(1,i)    = acos(aux_y/sqrt(aux_y**2.+aux_z**2.))
         !  angpos(2,i)    = acos(aux_z/sqrt(aux_x**2.+aux_z**2.))
         !  angpos(3,i)    = acos(aux_x/sqrt(aux_x**2.+aux_y**2.))! or ap(p)%phi
         !  angvel(1,i)    = ap(p)%omx 
         !  angvel(2,i)    = ap(p)%omy
         !  angvel(3,i)    = ap(p)%omz
          fltot(1,i)     = rkp(p)%fxltot
          fltot(2,i)     = rkp(p)%fyltot
          fltot(3,i)     = rkp(p)%fzltot
          accel(1,i)     = rkp(p)%dudt 
          accel(2,i)     = rkp(p)%dvdt
          accel(3,i)     = rkp(p)%dwdt
         !  angaccel(1,i)  = rkp(p)%domxdt 
         !  angaccel(2,i)  = rkp(p)%domydt
         !  angaccel(3,i)  = rkp(p)%domzdt
         !  extravec(1,i)  = radius*aux_x 
         !  extravec(2,i)  = radius*aux_y 
         !  extravec(3,i)  = radius*aux_z
         !  extrascal(1,i) = 0.
       endif
    enddo
    skipacc = 0
    call MPI_FILE_OPEN(MPI_COMM_WORLD, 'partvisu'//istepchar//'.bin', &
         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, error)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,error)  ! guarantee overwriting
    !skip = 3 ! vector
    !disp = np*skipacc*lenr + mydisp*skip*lenr
    !call MPI_FILE_SET_VIEW(fh, disp, MPI_REAL8,MPI_REAL8, 'native', & 
    !     MPI_INFO_NULL, error)
    !call MPI_FILE_WRITE(fh,posp(1,1),skip*npmstr,MPI_REAL8,MPI_STATUS_IGNORE,error)
    !skipacc = skipacc + skip
    !skip = 1 ! scalar
    !disp = np*skipacc*lenr + mydisp*skip*lenr
    !call MPI_FILE_SET_VIEW(fh, disp, MPI_REAL8,MPI_REAL8, 'native', & 
    !     MPI_INFO_NULL, error)
    !call MPI_FILE_WRITE(fh,ridp(1,1),skip*npmstr,MPI_REAL8,MPI_STATUS_IGNORE,error)
    !skipacc = skipacc + skip
    !skip = 3 ! vector
    !disp = np*skipacc*lenr + mydisp*skip*lenr
    !call MPI_FILE_SET_VIEW(fh, disp, MPI_REAL8,MPI_REAL8, 'native', & 
    !     MPI_INFO_NULL, error)
    !call MPI_FILE_WRITE(fh,velp(1,1),skip*npmstr,MPI_REAL8,MPI_STATUS_IGNORE,error)
    !skipacc = skipacc + skip
    !skip = 3 ! vector
    !disp = np*skipacc*lenr + mydisp*skip*lenr
    !call MPI_FILE_SET_VIEW(fh, disp, MPI_REAL8,MPI_REAL8, 'native', & 
    !     MPI_INFO_NULL, error)
    !call MPI_FILE_WRITE(fh,angpos(1,1),skip*npmstr,MPI_REAL8,MPI_STATUS_IGNORE,error)
    !skipacc = skipacc + skip
    ! Write force totals (fltot)
    skip = 3 ! vector
    disp = np*skipacc*lenr + mydisp*skip*lenr
    call MPI_FILE_SET_VIEW(fh, disp, MPI_REAL8,MPI_REAL8, 'native', & 
         MPI_INFO_NULL, error)
    call MPI_FILE_WRITE(fh,fltot(1,1),skip*npmstr,MPI_REAL8,MPI_STATUS_IGNORE,error)
    skipacc = skipacc + skip
    ! Write accelerations (accel)
    skip = 3 ! vector
    disp = np*skipacc*lenr + mydisp*skip*lenr
    call MPI_FILE_SET_VIEW(fh, disp, MPI_REAL8,MPI_REAL8, 'native', & 
         MPI_INFO_NULL, error)
    call MPI_FILE_WRITE(fh,accel(1,1),skip*npmstr,MPI_REAL8,MPI_STATUS_IGNORE,error)
    skipacc = skipacc + skip
    !skipacc = skipacc + skip
    !skip = 3 ! vector
    !disp = np*skipacc*lenr + mydisp*skip*lenr
    !call MPI_FILE_SET_VIEW(fh, disp, MPI_REAL8,MPI_REAL8, 'native', & 
    !     MPI_INFO_NULL, error)
    !call MPI_FILE_WRITE(fh,angaccel(1,1),skip*npmstr,MPI_REAL8,MPI_STATUS_IGNORE,error)
    !skipacc = skipacc + skip
    !skip = 3 ! vector
    !disp = np*skipacc*lenr + mydisp*skip*lenr
    !call MPI_FILE_SET_VIEW(fh, disp, MPI_REAL8,MPI_REAL8, 'native', & 
    !     MPI_INFO_NULL, error)
    !call MPI_FILE_WRITE(fh,extravec(1,1),skip*npmstr,MPI_REAL8,MPI_STATUS_IGNORE,error)
    !skipacc = skipacc + skip
    !skip = 1 ! scalar
    !disp = np*skipacc*lenr + mydisp*skip*lenr
    !call MPI_FILE_SET_VIEW(fh, disp, MPI_REAL8,MPI_REAL8, 'native', & 
    !     MPI_INFO_NULL, error)
    !call MPI_FILE_WRITE(fh,extrascal(1,1),skip*npmstr,MPI_REAL8,MPI_STATUS_IGNORE,error)
    call MPI_FILE_CLOSE(fh,error)
    !deallocate(posp)
    !deallocate(velp)
    !deallocate(ridp)
    !deallocate(angpos)
    !deallocate(angvel)
    deallocate(fltot)
    deallocate(accel)
    !deallocate(angaccel)
    !deallocate(extravec)
    !deallocate(extrascal)
    !
    return
  end subroutine outputpart
  !
  subroutine historypart(nr,icount)
    implicit none
    integer :: skip,skipacc
    integer,dimension(0:dims(1)*dims(2)-1) :: npmstr_glob,npmstr_glob_all
    real, allocatable, dimension(:,:) :: posp,velp,ridp,angpos,angvel,accel,angaccel,extravec,extrascal 
    real :: aux_x,aux_y,aux_z
    integer i,j,p,idp
    integer :: fh
    integer,intent(in) :: nr,icount
    real :: xdum
    integer :: lenr
    integer(kind=MPI_OFFSET_KIND) :: disp
    integer :: mydisp
    character(len=7) :: istepchar
    type save_history
       real :: idp,x,y,z,u,v,w
    end type save_history
    type(save_history) :: hp
    integer :: size_history
    !
    !inquire (iolength=lenr) xdum
    lenr = sizeof(xdum)
    size_history = 7
    write(istepchar,'(i7.7)') nr
    i = 0
    call MPI_FILE_OPEN(MPI_COMM_WORLD, 'parthist'//istepchar//'.bin', &
         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, error)
    do p=1,pmax
       if(ap(p)%mslv.gt.0) then
          hp%idp = 1.*ap(p)%mslv 
          hp%x = ap(p)%x 
          hp%y = ap(p)%y 
          hp%z = ap(p)%z 
          hp%u = ap(p)%u 
          hp%v = ap(p)%v 
          hp%w = ap(p)%w 
          idp  = ap(p)%mslv
          disp = ((idp-1)+(np-1)*(icount-1))*size_history*lenr
          call MPI_FILE_SET_VIEW(fh, disp, MPI_REAL8,MPI_REAL8, 'native', & 
               MPI_INFO_NULL, error)
          call MPI_FILE_WRITE(fh,hp%x,size_history,MPI_REAL8,MPI_STATUS_IGNORE,error)
       endif
    enddo
    call MPI_FILE_CLOSE(fh,error)
    return
  end subroutine historypart
  !
end module mod_output

