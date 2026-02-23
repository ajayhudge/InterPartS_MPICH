module mod_autocorr
  use decomp_2d
  use mod_param
  use mod_common
  use mod_common_mpi
  use mod_fftw
  use mod_phase_indicator
contains
  !
  subroutine autocorr(istep)
    implicit none
    integer, intent(in) :: istep
    integer :: ii,i,j,k,kglob
    real, allocatable, dimension(:,:,:) :: ux,vx,wx,uy,vy,wy,uz,vz,wz
    real, allocatable, dimension(:,:,:) :: gamu,gamv,gamw,gamp,upart,vpart,wpart
    real, allocatable, dimension(:,:) :: Ruuj,Rvvj,Rwwj,Rglbj,sum_allj
    real, allocatable, dimension(:,:) :: Ruui,Rvvi,Rwwi,Rglbi,sum_alli
    complex, allocatable, dimension(:) :: arrcmplx
    real,    allocatable, dimension(:) :: arrreal
    integer(kind=MPI_OFFSET_KIND) :: filesize,disp
    real :: xdum
    integer :: fh,lenr
    integer :: counter,group_world,group_const_z_ypencil,comm_const_z_ypencil, &
         group_const_z_xpencil,comm_const_z_xpencil
    integer, dimension(dims(1)) :: ranks_const_z_ypencil,ranks_const_z_xpencil
    character(len=7) :: istepchar
    real, dimension(ktot) :: um,vm,wm,sum_all
    !
    write(istepchar,'(i7.7)') istep
    !
    ! spanwise    direction -> i
    ! streamwise  direction -> j
    ! wall-normal direction -> k
    !
    lenr = sizeof(xdum)
    !
    ! initialize groups of processes with the same
    ! wall-normal coordinate for averaging purposes
    !
    ! y-pencils
    !
    call MPI_COMM_GROUP(MPI_COMM_WORLD,group_world,error)
    counter = 0
    do i=0,dims(1)-1
       do k=0,dims(2)-1
          if((ystart(3)-1)/(ktot/dims(2)).eq.k) then
             ranks_const_z_ypencil(i+1) = counter
          endif
          counter = counter + 1
       enddo
    enddo
    call MPI_GROUP_INCL(group_world, dims(1), ranks_const_z_ypencil, group_const_z_ypencil, error)
    call MPI_COMM_CREATE(MPI_COMM_WORLD, group_const_z_ypencil, comm_const_z_ypencil,error)  
    !
    ! x-pencils
    !
    counter = 0
    do j=0,dims(1)-1
       do k=0,dims(2)-1
          if((xstart(3)-1)/(ktot/dims(2)).eq.k) then
             ranks_const_z_xpencil(j+1) = counter
          endif
          counter = counter + 1
       enddo
    enddo
    call MPI_GROUP_INCL(group_world, dims(1), ranks_const_z_xpencil, group_const_z_xpencil, error)
    call MPI_COMM_CREATE(MPI_COMM_WORLD, group_const_z_xpencil, comm_const_z_xpencil,error)  
    allocate(gamu(0:i1,0:j1,0:k1), &
         gamv(0:i1,0:j1,0:k1), &
         gamw(0:i1,0:j1,0:k1), &
         gamp(0:i1,0:j1,0:k1), &
         upart(0:i1,0:j1,0:k1), &
         vpart(0:i1,0:j1,0:k1), &
         wpart(0:i1,0:j1,0:k1) )
    call phase_indicator(gamu,gamv,gamw,gamp,upart,vpart,wpart,1)
    !
    ! impose rigid body motion in the flow-field
    !
    allocate(uz(imax,jmax,kmax))
    allocate(vz(imax,jmax,kmax))
    allocate(wz(imax,jmax,kmax))
    !$omp workshare
    uz(:,:,:) = unew(1:imax,1:jmax,1:kmax)*gamu(1:imax,1:jmax,1:kmax)+upart(1:imax,1:jmax,1:kmax)*(1.-gamu(1:imax,1:jmax,1:kmax))
    vz(:,:,:) = vnew(1:imax,1:jmax,1:kmax)*gamv(1:imax,1:jmax,1:kmax)+vpart(1:imax,1:jmax,1:kmax)*(1.-gamv(1:imax,1:jmax,1:kmax))
    wz(:,:,:) = wnew(1:imax,1:jmax,1:kmax)*gamw(1:imax,1:jmax,1:kmax)+wpart(1:imax,1:jmax,1:kmax)*(1.-gamw(1:imax,1:jmax,1:kmax))
    !$omp end workshare
!
    !$omp parallel default(none) &
    !$omp&shared(uz,vz,wz) &
    !$omp&private(i,j,k) &
    !$omp&reduction(+:um,vm,wm)
    !$omp do
    do k=1,kmax
       um(k) = 0.
       vm(k) = 0.
       wm(k) = 0.
       do j=1,jmax
          do i=1,imax
             um(k)  = um(k)  + uz(i,j,k)
             vm(k)  = vm(k)  + vz(i,j,k)
             wm(k)  = wm(k)  + wz(i,j,k)
          enddo
       enddo
    enddo
    !$omp end parallel
    !
    call mpi_allreduce(um(1),sum_all(1),ktot,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    um(:) = sum_all(:)/(1.*itot*jtot)
    !$omp end workshare
    call mpi_allreduce(vm(1),sum_all(1),ktot,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    vm(:) = sum_all(:)/(1.*itot*jtot)
    !$omp end workshare
    call mpi_allreduce(wm(1),sum_all(1),ktot,mpi_real8,mpi_sum,comm_cart,error)
    !$omp workshare
    wm(:) = sum_all(:)/(1.*itot*jtot)
    !$omp end workshare
    !$omp parallel default(none) &
    !$omp&shared(uz,vz,wz,um,vm,wm) &
    !$omp&private(i,j,k)
    !$omp do
    do k=1,kmax
       do j=1,jmax
          do i=1,imax
             uz(i,j,k) = uz(i,j,k) - um(k)
             vz(i,j,k) = vz(i,j,k) - vm(k)
             wz(i,j,k) = wz(i,j,k) - wm(k)
          enddo
       enddo
    enddo
    !$omp end parallel
    !
    deallocate(gamu,gamv,gamw,gamp,upart,vpart,wpart)
    !
    allocate(uy(itot/dims(1),jtot,ktot/dims(2)))
    allocate(vy(itot/dims(1),jtot,ktot/dims(2)))
    allocate(wy(itot/dims(1),jtot,ktot/dims(2)))
    !
    ! align computational topology in the streamwise direction
    ! (y-pencils)
    !
    call transpose_z_to_y(uz,uy)
    call transpose_z_to_y(vz,vy)
    call transpose_z_to_y(wz,wy)
    deallocate(uz,vz,wz)
    !
    ! compute autocorrelation for streamwise
    ! separation distances
    !
    allocate(arrcmplx(jtot/2+1))
    allocate(arrreal(jtot))
    allocate(Ruuj(jtot,ktot/dims(2)))
    allocate(Rvvj(jtot,ktot/dims(2)))
    allocate(Rwwj(jtot,ktot/dims(2)))
    !$omp workshare
    arrcmplx(:) = 0.
    arrreal(:)  = 0.
    !$omp end workshare
    !  call init_fft(jtot)
    !$omp parallel default(none) &
    !$omp&shared(uy,vy,wy,plan_r2c_y,plan_c2r_y,ysize) &
    !$omp&private(i,j,k,kglob,arrreal,arrcmplx) &
    !$omp&reduction(+:Ruuj,Rvvj,Rwwj) 
    !$omp do
    do k=1,ktot/dims(2)
       Ruuj(:,k) = 0.
       Rvvj(:,k) = 0.
       Rwwj(:,k) = 0.
       do i=1,itot/dims(1)
          kglob = k + ysize(3)-1
          call fftr2c(jtot,uy(i,:,k),arrcmplx(:),plan_r2c_y)
          arrcmplx(:) = arrcmplx(:)*conjg(arrcmplx(:))
          call fftc2r(jtot,arrcmplx(:),arrreal(:),plan_c2r_y)
          Ruuj(:,k) = Ruuj(:,k) + arrreal(:)
          !
          call fftr2c(jtot,vy(i,:,k),arrcmplx(:),plan_r2c_y)
          arrcmplx(:) = arrcmplx(:)*conjg(arrcmplx(:))
          call fftc2r(jtot,arrcmplx(:),arrreal(:),plan_c2r_y)
          Rvvj(:,k) = Rvvj(:,k) + arrreal(:)
          !
          call fftr2c(jtot,wy(i,:,k),arrcmplx(:),plan_r2c_y)
          arrcmplx(:) = arrcmplx(:)*conjg(arrcmplx(:))
          call fftc2r(jtot,arrcmplx(:),arrreal(:),plan_c2r_y)
          Rwwj(:,k) = Rwwj(:,k) + arrreal(:)
       enddo
    enddo
    !$omp end parallel
    deallocate(arrreal,arrcmplx)
    !  call clean_fft
    !
    allocate(sum_allj(jtot,ktot/dims(2)))
    call mpi_allreduce(Ruuj(1,1),sum_allj(1,1),jtot*ktot/dims(2),MPI_REAL8,MPI_SUM,comm_const_z_ypencil,error)
    !$omp workshare
    Ruuj(:,:) = sum_allj(:,:)/(1.*itot)
    !$omp end workshare
    call mpi_allreduce(Rvvj(1,1),sum_allj(1,1),jtot*ktot/dims(2),MPI_REAL8,MPI_SUM,comm_const_z_ypencil,error)
    !$omp workshare
    Rvvj(:,:) = sum_allj(:,:)/(1.*itot)
    !$omp end workshare
    call mpi_allreduce(Rwwj(1,1),sum_allj(1,1),jtot*ktot/dims(2),MPI_REAL8,MPI_SUM,comm_const_z_ypencil,error)
    !$omp workshare
    Rwwj(:,:) = sum_allj(:,:)/(1.*itot)
    !$omp end workshare
    !
    allocate(Rglbj(jtot/2,ktot/dims(2)))
    !
    ! prepare and write data with mpi-i/o
    !
    !$omp parallel default(none) &
    !$omp&shared(Rglbj,Ruuj) &
    !$omp&private(j,k)
    !$omp do
    do k=1,ktot/dims(2)
       Rglbj(1,k) = Ruuj(1,k)/Ruuj(1,k)
       if(Rglbj(1,k).ne.Rglbj(1,k)) Rglbj(1,k) = 1. ! to prevent NaN at the wall
       do j = 1,jtot/2-1
          Rglbj(j+1,k) = 0.5*(Ruuj(j+1,k)+Ruuj(jtot-j+1,k))/Ruuj(1,k)
          if(Rglbj(j+1,k).ne.Rglbj(j+1,k)) Rglbj(1,k) = 1. ! to prevent NaN at the wall
       enddo
    enddo
    !$omp end parallel
    call MPI_FILE_OPEN(MPI_COMM_WORLD, 'ruuj_fld_'//istepchar//'_2d.bin', &
         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, error)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,error)  ! guarantee overwriting
    disp = ((ystart(3)-1)/(ktot/dims(2)))*jtot/2*ktot/dims(2)*lenr
    call MPI_FILE_SET_VIEW(fh, disp, MPI_REAL8,MPI_REAL8, 'native', &
         MPI_INFO_NULL, error)
    call MPI_FILE_WRITE(fh,Rglbj(1,1),jtot/2*ktot/dims(2),MPI_REAL8,MPI_STATUS_IGNORE,error)
    call MPI_FILE_CLOSE(fh,error)
    !
    !$omp parallel default(none) &
    !$omp&shared(Rglbj,Rvvj) &
    !$omp&private(j,k)
    !$omp do
    do k=1,ktot/dims(2)
       Rglbj(1,k) = Rvvj(1,k)/Rvvj(1,k)
       if(Rglbj(1,k).ne.Rglbj(1,k)) Rglbj(1,k) = 1. ! to prevent NaN at the wall
       do j = 1,jtot/2-1
          Rglbj(j+1,k) = 0.5*(Rvvj(j+1,k)+Rvvj(jtot-j+1,k))/Rvvj(1,k)
          if(Rglbj(j+1,k).ne.Rglbj(j+1,k)) Rglbj(1,k) = 1. ! to prevent NaN at the wall
       enddo
    enddo
    !$omp end parallel
    call MPI_FILE_OPEN(MPI_COMM_WORLD, 'rvvj_fld_'//istepchar//'_2d.bin', &
         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, error)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,error)  ! guarantee overwriting
    disp = ((ystart(3)-1)/(ktot/dims(2)))*jtot/2*ktot/dims(2)*lenr
    call MPI_FILE_SET_VIEW(fh, disp, MPI_REAL8,MPI_REAL8, 'native', &
         MPI_INFO_NULL, error)
    call MPI_FILE_WRITE(fh,Rglbj(1,1),jtot/2*ktot/dims(2),MPI_REAL8,MPI_STATUS_IGNORE,error)
    call MPI_FILE_CLOSE(fh,error)
    !
    !$omp parallel default(none) &
    !$omp&shared(Rglbj,Rwwj) &
    !$omp&private(j,k)
    !$omp do
    do k=1,ktot/dims(2)
       Rglbj(1,k) = Rwwj(1,k)/Rwwj(1,k)
       if(Rglbj(1,k).ne.Rglbj(1,k)) Rglbj(1,k) = 1. ! to prevent NaN at the wall
       do j = 1,jtot/2-1
          Rglbj(j+1,k) = 0.5*(Rwwj(j+1,k)+Rwwj(jtot-j+1,k))/Rwwj(1,k)
          if(Rglbj(j+1,k).ne.Rglbj(j+1,k)) Rglbj(j+1,k) = 1. ! to prevent NaN at the wall
       enddo
    enddo
    !$omp end parallel
    call MPI_FILE_OPEN(MPI_COMM_WORLD, 'rwwj_fld_'//istepchar//'_2d.bin', &
         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, error)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,error)  ! guarantee overwriting
    disp = ((ystart(3)-1)/(ktot/dims(2)))*jtot/2*ktot/dims(2)*lenr
    call MPI_FILE_SET_VIEW(fh, disp, MPI_REAL8,MPI_REAL8, 'native', &
         MPI_INFO_NULL, error)
    call MPI_FILE_WRITE(fh,Rglbj(1,1),jtot/2*ktot/dims(2),MPI_REAL8,MPI_STATUS_IGNORE,error)
    call MPI_FILE_CLOSE(fh,error)
    deallocate(Rglbj,sum_allj)
    deallocate(Ruuj,Rvvj,Rwwj)
    !
    ! autocorrelation in the streamwise direction done
    !
    !
    ! align computational topology in the spanwise direction
    ! (x-pencils)
    !
    allocate(ux(itot,jtot/dims(1),1:ktot/dims(2)))
    allocate(vx(itot,jtot/dims(1),1:ktot/dims(2)))
    allocate(wx(itot,jtot/dims(1),1:ktot/dims(2)))
    call transpose_y_to_x(uy,ux)
    call transpose_y_to_x(vy,vx)
    call transpose_y_to_x(wy,wx)
    deallocate(uy,vy,wy)
    !
    ! compute autocorrelation for spanwise
    ! separation distances
    !
    allocate(arrcmplx(itot/2+1))
    allocate(arrreal(itot))
    allocate(Ruui(itot,ktot/dims(2)))
    allocate(Rvvi(itot,ktot/dims(2)))
    allocate(Rwwi(itot,ktot/dims(2)))
    !  call init_fft(itot)
    !$omp parallel default(none) &
    !$omp&shared(ux,vx,wx,plan_r2c_x,plan_c2r_x,xsize) &
    !$omp&private(i,j,k,kglob,arrreal,arrcmplx) &
    !$omp&reduction(+:Ruui,Rvvi,Rwwi) 
    !$omp do
    do k=1, ktot/dims(2)
       Ruui(:,k) = 0.
       Rvvi(:,k) = 0.
       Rwwi(:,k) = 0.
       do j=1, jtot/dims(1)
          kglob = k + xsize(3)-1
          call fftr2c(itot,ux(:,j,k),arrcmplx(:),plan_r2c_x)
          arrcmplx(:) = arrcmplx(:)*conjg(arrcmplx(:))
          call fftc2r(itot,arrcmplx(:),arrreal(:),plan_c2r_x)
          Ruui(:,k) = Ruui(:,k) + arrreal(:)
          !
          call fftr2c(itot,vx(:,j,k),arrcmplx(:),plan_r2c_x)
          arrcmplx(:) = arrcmplx(:)*conjg(arrcmplx(:))
          call fftc2r(itot,arrcmplx(:),arrreal(:),plan_c2r_x)
          Rvvi(:,k) = Rvvi(:,k) + arrreal(:)
          !
          call fftr2c(itot,wx(:,j,k),arrcmplx(:),plan_r2c_x)
          arrcmplx(:) = arrcmplx(:)*conjg(arrcmplx(:))
          call fftc2r(itot,arrcmplx(:),arrreal(:),plan_c2r_x)
          Rwwi(:,k) = Rwwi(:,k) + arrreal(:)
       enddo
    enddo
    !$omp end parallel
    deallocate(ux,vx,wx)
    deallocate(arrreal,arrcmplx)
    !  call clean_fft
    !
    allocate(sum_alli(itot,ktot/dims(2)))
    call mpi_allreduce(Ruui(1,1),sum_alli(1,1),itot*ktot/dims(2),MPI_REAL8,MPI_SUM,comm_const_z_xpencil,error)
    !$omp workshare
    Ruui(:,:) = sum_alli(:,:)/(1.*jtot)
    !$omp end workshare
    call mpi_allreduce(Rvvi(1,1),sum_alli(1,1),itot*ktot/dims(2),MPI_REAL8,MPI_SUM,comm_const_z_xpencil,error)
    !$omp workshare
    Rvvi(:,:) = sum_alli(:,:)/(1.*jtot)
    !$omp end workshare
    call mpi_allreduce(Rwwi(1,1),sum_alli(1,1),itot*ktot/dims(2),MPI_REAL8,MPI_SUM,comm_const_z_xpencil,error)
    !$omp workshare
    Rwwi(:,:) = sum_alli(:,:)/(1.*jtot)
    !$omp end workshare
    !
    ! prepare and write data with mpi-i/o
    !
    allocate(Rglbi(itot/2,ktot/dims(2)))
    !
    !$omp parallel default(none) &
    !$omp&shared(Rglbi,Ruui) &
    !$omp&private(i,k)
    !$omp do
    do k=1,ktot/dims(2)
       Rglbi(1,k) = Ruui(1,k)/Ruui(1,k)
       if(Rglbi(1,k).ne.Rglbi(1,k)) Rglbi(1,k) = 1. ! to prevent NaN at the wall
       do i = 1,itot/2-1
          Rglbi(i+1,k) = 0.5*(Ruui(i+1,k)+Ruui(itot-i+1,k))/Ruui(1,k)
          if(Rglbi(i+1,k).ne.Rglbi(i+1,k)) Rglbi(1,k) = 1. ! to prevent NaN at the wall
       enddo
    enddo
    !$omp end parallel
    call MPI_FILE_OPEN(MPI_COMM_WORLD, 'ruui_fld_'//istepchar//'_2d.bin', &
         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, error)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,error)  ! guarantee overwriting
    disp = ((xstart(3)-1)/(ktot/dims(2)))*itot/2*ktot/dims(2)*lenr
    call MPI_FILE_SET_VIEW(fh, disp, MPI_REAL8,MPI_REAL8, 'native', &
         MPI_INFO_NULL, error)
    call MPI_FILE_WRITE(fh,Rglbi(1,1),itot/2*ktot/dims(2),MPI_REAL8,MPI_STATUS_IGNORE,error)
    call MPI_FILE_CLOSE(fh,error)
    !
    !$omp parallel default(none) &
    !$omp&shared(Rglbi,Rvvi) &
    !$omp&private(i,k)
    !$omp do
    do k=1,ktot/dims(2)
       Rglbi(1,k) = Rvvi(1,k)/Rvvi(1,k)
       if(Rglbi(1,k).ne.Rglbi(1,k)) Rglbi(1,k) = 1. ! to prevent NaN at the wall
       do i = 1,itot/2-1
          Rglbi(i+1,k) = 0.5*(Rvvi(i+1,k)+Rvvi(itot-i+1,k))/Rvvi(1,k)
          if(Rglbi(i+1,k).ne.Rglbi(i+1,k)) Rglbi(1,k) = 1. ! to prevent NaN at the wall
       enddo
    enddo
    !$omp end parallel
    call MPI_FILE_OPEN(MPI_COMM_WORLD, 'rvvi_fld_'//istepchar//'_2d.bin', &
         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, error)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,error)  ! guarantee overwriting
    disp = ((xstart(3)-1)/(ktot/dims(2)))*itot/2*ktot/dims(2)*lenr
    call MPI_FILE_SET_VIEW(fh, disp, MPI_REAL8,MPI_REAL8, 'native', &
         MPI_INFO_NULL, error)
    call MPI_FILE_WRITE(fh,Rglbi(1,1),itot/2*ktot/dims(2),MPI_REAL8,MPI_STATUS_IGNORE,error)
    call MPI_FILE_CLOSE(fh,error)
    !
    !$omp parallel default(none) &
    !$omp&shared(Rglbi,Rwwi) &
    !$omp&private(i,k)
    !$omp do
    do k=1,ktot/dims(2)
       Rglbi(1,k) = Rwwi(1,k)/Rwwi(1,k)
       if(Rglbi(1,k).ne.Rglbi(1,k)) Rglbi(1,k) = 1. ! to prevent NaN at the wall
       do i = 1,itot/2-1
          Rglbi(i+1,k) = 0.5*(Rwwi(i+1,k)+Rwwi(itot-i+1,k))/Rwwi(1,k)
          if(Rglbi(i+1,k).ne.Rglbi(i+1,k)) Rglbi(i+1,k) = 1. ! to prevent NaN at the wall
       enddo
    enddo
    !$omp end parallel
    call MPI_FILE_OPEN(MPI_COMM_WORLD, 'rwwi_fld_'//istepchar//'_2d.bin', &
         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, error)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,error)  ! guarantee overwriting
    disp = ((xstart(3)-1)/(ktot/dims(2)))*itot/2*ktot/dims(2)*lenr
    call MPI_FILE_SET_VIEW(fh, disp, MPI_REAL8,MPI_REAL8, 'native', &
         MPI_INFO_NULL, error)
    call MPI_FILE_WRITE(fh,Rglbi(1,1),itot/2*ktot/dims(2),MPI_REAL8,MPI_STATUS_IGNORE,error)
    call MPI_FILE_CLOSE(fh,error)
    deallocate(Rglbi,sum_alli)
    deallocate(Ruui,Rvvi,Rwwi)
    !
    ! autocorrelation in the spanwise direction done
    !
    return
  end subroutine autocorr
  !
end module mod_autocorr
