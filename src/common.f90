module mod_common
  use mod_param
  !
  real ,dimension(0:i1,0:j1,0:k1) :: unew,vnew,wnew,pnew,cnew, &
       dudtold,dvdtold,dwdtold, dcdtold
  real, dimension(-1:i1+1,-1:j1+1,-1:k1+1) :: dudtf,dvdtf,dwdtf, dcdtf
  real, dimension(0:i1,0:j1,0:k1) :: dudt,dvdt,dwdt, dcdt, source 
  real, dimension(-1:i1+1,-1:j1+1,-1:k1+1) :: forcex,forcey,forcez
  real(mytype) :: time,dt
  real(mytype) ::  wi(itot+15), wj(jtot+15)
  real, dimension(imax,jmax) :: xyrt
  real, dimension(kmax) :: a,b,c
  real :: forcextot,forceytot,forceztot
  real :: u_bulk,v_bulk,w_bulk
  real wallshearold,wallshearnew
  real :: dpdx_sumrk,dpdy_new
  integer :: rkiter
  real :: rkparalpha
  !
  ! particles
  !
  ! for simplifying the communication between threads when there is a
  ! new master, the derived type 'particle' should be organized in this way:
  !
  !   (1st) real data that has to be communicated when there is a new master
  !         is defined contiguously
  !
  !   (2nd) integer data that has to be communicated when there is a new master
  !         is defined contiguously
  !
  !   (3rd) real data that does not have to be communicated
  !
  !   (4th) integer data that does not have to be communicated
  !
  type particle
     real :: x,y,z,theta,phi, &
          u,v,w, &
          omx,omy,omz,omtheta, &
          intu,intv,intw, &
          intomx,intomy,intomz, &
          colfx,colfy,colfz, &
          coltx,colty,coltz ! 24
     real, dimension(nqmax) :: dx,dy,dz, &
          dxt,dyt,dzt, & 
          dut,dvt,dwt, &
          psi ! 10*nqmax
     ! total ammount of reals to be communicated: 24+10*nqmax
     integer :: qmax ! 1
     integer, dimension(nqmax) :: firstc ! 1*nqmax
     ! total ammount of integers to be communicated: 1+nqmax
     integer, dimension(nfriendsmax) :: friend ! 1*nfriends
     integer :: nfriends
     real :: fxltot,fyltot,fzltot, &
          torqxltot,torqyltot,torqzltot,torqtheta
     real, dimension(nlmax) :: xfp,yfp,zfp, &
          ul,vl,wl, &
          dudtl,dvdtl,dwdtl, &
          fxl,fyl,fzl
     real :: vol,mominert,ratiorho
     integer :: mslv ! master/slave flag 
     ! > 0: master, < 0: slave, = 0: not active (SS)
     ! Master processes collect the integral contributions of all slave processes. (SS)
     ! Slave processes compute partial integrals within their own domains. (SS)
     ! Distributed spherical integration is implemented via MPI communication. (SS)
     ! integer, dimension(8) :: nb !! dangerous accesses at nb(0) 
     integer, dimension(0:8) :: nb
     ! index of the neighbor to which the particle has to be sent (SS)
     !  6---7---8
     !  5---0---1  (0 is the current process, 1-8 are neighbors)
     !  4---3---2  (SS)
  end type particle
  type(particle), dimension(npmax) :: ap ! 'a particle' array
  ! this contains all the particle info
  !
  type particle_old
     real :: x,y,z,theta,phi, &
          u,v,w, &
          omx,omy,omz,omtheta, &
          intu,intv,intw, &
          intomx,intomy,intomz, &
          colfx,colfy,colfz, &
          coltx,colty,coltz ! 24
     real, dimension(nqmax) :: dx,dy,dz, &
          dxt,dyt,dzt, &
          dut,dvt,dwt, &
          psi ! 10*nqmax
     integer :: qmax ! 1
     integer, dimension(nqmax) :: firstc ! 1*nqmax
     integer :: mslv
  end type particle_old
  type(particle_old), dimension(npmax) :: op !old particle array (integration of N-E equations)
  type(particle_old), dimension(npmax) :: sp !send particle array (re-ordering of masters)
  !
  ! The structure of 'particle_sumrk' follows the same criterion as 'particle'
  ! for its organization
  !
  type particle_sumrk
     real :: fxltot,fyltot,fzltot, &
             torqxltot,torqyltot,torqzltot,torqtheta, &
             colfx,colfy,colfz, &
             coltx,colty,coltz, &
             dudt,dvdt,dwdt, &
             domxdt,domydt,domzdt
  end type particle_sumrk
  type(particle_sumrk), dimension(npmax) :: rkp,srkp ! this contains some particle
  ! data integrated over the rk3 substeps
  type particle_interior
     real :: x,y,z
     !real, dimension(nl+nl2+nl3+nl4) :: xfp,yfp,zfp,dvlagr
     real, dimension(nltotmax) :: xfp,yfp,zfp,dvlagr
  end type particle_interior
  !
  real :: dVlagr,dVeul
  real :: radfp
  real, dimension(NLmax) :: thetarc,phirc ! angular position of the Lfps
  !
  integer :: pmax,npmstr
  integer, dimension(npmax) :: nla
  !
end module mod_common
!
module mod_common_mpi
  use mpi
  use decomp_2d
  implicit none
  integer :: myid,xhalo,yhalo,restartp,rankcw
  integer :: comm_cart
  integer :: xhalo_ibm,yhalo_ibm
  !
  logical periods(3),reorder
  integer error,request,status(MPI_STATUS_SIZE)
  integer right,rightfront,front,leftfront,left,leftback,back,rightback
  integer, dimension(0:8) :: neighbor
  integer, dimension(1:2) :: coords
  real :: boundleftmyid,boundfrontmyid
  real :: presgrad
  !
end module mod_common_mpi
