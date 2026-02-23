module mod_coordsfp
  use mod_common
  use mod_common_mpi
  use mod_param
  implicit none
  private
  public coordsfp,coordsfp_interior
contains
  !
  subroutine coordsfp
    implicit none
    integer :: l,p,ll
    real :: dummy
    real :: xfploc,yfploc,zfploc
    real :: isperiodx,isperiody
    real :: coorxfp,cooryfp,coorzfp
    logical :: isout
    !
    ! position of Lfps wrt the center of sphere
    !
    open(25,file='lagrangianforcepoints2')
    !open(25,file='lagrangianforcepoints2.inp')
    read(25,*)
    read(25,*)
    read(25,*)
    do l=1,NL
       read(25,'(6E16.8)') dummy,dummy,dummy,thetarc(l),phirc(l),radfp
    enddo
    close(25)
    !
    ! position of Lfps wrt to the center of computational domain
    !
    !$omp workshare
    nla(:) = 0 ! set to zero from 1 to npmax
    !$omp end workshare
!
    !$omp parallel default(none)               &
    !$omp&shared(ap,nla,pmax)                  &
    !$omp&shared(radfp,phirc,thetarc)          &
    !$omp&shared(boundleftmyid,boundfrontmyid) &
    !$omp&private(p,l,ll,coorxfp,cooryfp,xfploc,yfploc,zfploc)    &
    !$omp&private(isperiodx,isperiody,isout)
    !$omp do 
    do p=1,pmax
       if (ap(p)%mslv .ne. 0) then
          ! myid is master of particle ap(p)%mslv
          ll = 0
          do l=1,NL
             xfploc = ap(p)%x + radfp*sin(ap(p)%theta*0.+thetarc(l))*cos(ap(p)%phi*0.+phirc(l))
             yfploc = ap(p)%y + radfp*sin(ap(p)%theta*0.+thetarc(l))*sin(ap(p)%phi*0.+phirc(l))
             zfploc = ap(p)%z + radfp*cos(ap(p)%theta*0.+thetarc(l))
             isperiodx = 0.
             isperiody = 0.
             if (xfploc.lt.0.+0.5*dx) isperiodx =  1.
             if (xfploc.ge.lx+0.5*dx) isperiodx = -1.
             if (yfploc.lt.0.+0.5*dy) isperiody =  1.
             if (yfploc.ge.ly+0.5*dy) isperiody = -1.
             isout = .false.
             coorxfp = (xfploc+isperiodx*lx-boundleftmyid )*dxi
             if( nint(coorxfp).lt.1 .or. nint(coorxfp) .gt.imax ) isout = .true.
             cooryfp = (yfploc+isperiody*ly-boundfrontmyid)*dyi
             if( nint(cooryfp).lt.1 .or. nint(cooryfp) .gt.jmax ) isout = .true.
!
             if(.not.isout) then
               ll = ll + 1
               ap(p)%xfp(ll) = xfploc
               ap(p)%yfp(ll) = yfploc
               ap(p)%zfp(ll) = zfploc
               ap(p)%ul(ll) = ap(p)%u + ap(p)%omy*(ap(p)%zfp(ll)-ap(p)%z) &
                                      - ap(p)%omz*(ap(p)%yfp(ll)-ap(p)%y)
               ap(p)%vl(ll) = ap(p)%v + ap(p)%omz*(ap(p)%xfp(ll)-ap(p)%x) &
                                      - ap(p)%omx*(ap(p)%zfp(ll)-ap(p)%z)
               ap(p)%wl(ll) = ap(p)%w + ap(p)%omx*(ap(p)%yfp(ll)-ap(p)%y) &
                                      - ap(p)%omy*(ap(p)%xfp(ll)-ap(p)%x)
             endif
          enddo
          nla(p) = ll
       endif
    enddo
    !$omp end parallel
    !
    ! dVlagr of outer shell
    !
    dVlagr = ( (4./3.)*pi*(radfp+0.5/dxi)**3. - &
         (4./3.)*pi*(radfp-0.5/dxi)**3. )/(1.*nl)
    !
    ! Volume of an Eulerian grid cell
    !
    dVeul  = 1./(dxi*dyi*dzi)
    !
    if (myid .eq. 0) then
       write(6,*) 'dVlagr, dVeul = ',dVlagr, dVeul
    endif
    !
    ! Define geometrical properties of the particles
    !
    !$omp parallel default(none)  &
    !$omp&shared(ap,nla,pmax)     &
    !$omp&private(p)
    !$omp do 
    do p=1,npmax
       ap(p)%vol = volp
       ap(p)%ratiorho = ratiorho
       ap(p)%mominert = mominert
    enddo
    !$omp end parallel
    !
    ! data to file
    !
    !if (ap(1)%mslv .eq. 1) then
    !  open(25,file=datadir//'lagrangianforcepoints')
    !  write(25,*) 'VARIABLES = "lfpx","lfpy","lfpz","radius"'
    !  write(25,*) 'ZONE T="Zone1"',' I=',NL,', F=POINT'
    !  write(25,*) ''
    !  do l=1,nl
    !    write(25,'(4E16.8)') ap(1)%xfp(l),ap(1)%yfp(l),ap(1)%zfp(l), &
    !                         sqrt( (ap(1)%xfp(l)-ap(1)%x)**2. + (ap(1)%yfp(l)-ap(1)%y)**2. + &
    !                               (ap(1)%zfp(l)-ap(1)%z)**2. )
    !  enddo
    !  close(25)
    !endif
    !
    return
  end subroutine coordsfp
  !
  subroutine coordsfp_interior(api)
    implicit none
    integer :: l,p,q
    type(particle_interior), dimension(npmax), intent(out) :: api
    real dummy,thetarcint(1:NLtot),phircint(1:NLtot),radfp2,radfp3,radfp4
    !
    ! position of Lagrangian force points wrt center of sphere
    !
    !open(42,file='poslfp/4ringstogether/data/lagrangianforcepoints2')
    open(42,file='lagrangianforcepoints2.inp')
    read(42,*)
    read(42,*)
    read(42,*)
    do l=1,NL
       read(42,'(6E16.8)') dummy,dummy,dummy,thetarc(l),phirc(l),radfp
       thetarcint(l) = thetarc(l)
       phircint(l)   = phirc(l)
    enddo
    do q=1,NL2
       l=NL+q
       read(42,'(6E16.8)') dummy,dummy,dummy,thetarcint(l),phircint(l),radfp2
    enddo
    do q=1,NL3
       l=NL+NL2+q
       read(42,'(6E16.8)') dummy,dummy,dummy,thetarcint(l),phircint(l),radfp3
    enddo
    do q=1,NL4
       l=NL+NL2+NL3+q
       read(42,'(6E16.8)') dummy,dummy,dummy,thetarcint(l),phircint(l),radfp4
    enddo
    close(42)
    !
    ! position of lfp's wrt to center of computational domain
    !
    do p=1,pmax
       if (ap(p)%mslv .gt. 0) then
          ! myid is master of particle ap(p)%mslv
          api(p)%x = ap(p)%x
          api(p)%y = ap(p)%y
          api(p)%z = ap(p)%z
          do l=1,NL
             api(p)%xfp(l)      = ap(p)%x + radfp*sin(ap(p)%theta+thetarc(l))*cos(ap(p)%phi+phirc(l))
             api(p)%yfp(l)      = ap(p)%y + radfp*sin(ap(p)%theta+thetarc(l))*sin(ap(p)%phi+phirc(l))
             api(p)%zfp(l)      = ap(p)%z + radfp*cos(ap(p)%theta+thetarc(l))
          enddo
          do q=1,NL2
             l=NL+q
             api(p)%xfp(l) = ap(p)%x + radfp2*sin(ap(p)%theta+thetarcint(l))*cos(ap(p)%phi+phircint(l))
             api(p)%yfp(l) = ap(p)%y + radfp2*sin(ap(p)%theta+thetarcint(l))*sin(ap(p)%phi+phircint(l))
             api(p)%zfp(l) = ap(p)%z + radfp2*cos(ap(p)%theta+thetarcint(l))
          enddo
          do q=1,NL3
             l=NL+NL2+q
             api(p)%xfp(l) = ap(p)%x + radfp3*sin(ap(p)%theta+thetarcint(l))*cos(ap(p)%phi+phircint(l))
             api(p)%yfp(l) = ap(p)%y + radfp3*sin(ap(p)%theta+thetarcint(l))*sin(ap(p)%phi+phircint(l))
             api(p)%zfp(l) = ap(p)%z + radfp3*cos(ap(p)%theta+thetarcint(l))
          enddo
          do q=1,NL4
             l=NL+NL2+NL3+q
             api(p)%xfp(l) = ap(p)%x + radfp4*sin(ap(p)%theta+thetarcint(l))*cos(ap(p)%phi+phircint(l))
             api(p)%yfp(l) = ap(p)%y + radfp4*sin(ap(p)%theta+thetarcint(l))*sin(ap(p)%phi+phircint(l))
             api(p)%zfp(l) = ap(p)%z + radfp4*cos(ap(p)%theta+thetarcint(l))
          enddo
       else
          ! myid is either slave of particle ap(p)%mslv or does not contain this particle
          do l=1,NL+nl2+nl3+nl4
             api(p)%xfp(l) = 0.
             api(p)%yfp(l) = 0.
             api(p)%zfp(l) = 0.
          enddo
       endif
    enddo
    !
    ! Volume of a Lagr. force point cell (lfp-cell)
    !
    do p=1,pmax
       do l=1,NL
          api(p)%dvlagr(l)= ( (4./3.)*pi*(radfp+0.5/dxi)**3. - &
               (4./3.)*pi*(radfp-0.5/dxi)**3. )/(1.*NL)
       enddo
       do q=1,NL2
          l=NL+q
          api(p)%dvlagr(l) = ( (4./3.)*pi*(radfp-0.5/dxi)**3. - &
               (4./3.)*pi*(radfp-1.5/dxi)**3. )/(1.*NL2)
       enddo
       do q=1,NL3
          l=NL+NL2+q
          api(p)%dvlagr(l) = ( (4./3.)*pi*(radfp-1.5/dxi)**3. - &
               (4./3.)*pi*(radfp-2.5/dxi)**3. )/(1.*NL3)
       enddo
       do q=1,NL4
          l=NL+NL2+NL3+q
          api(p)%dvlagr(l) = ( (4./3.)*pi*(radfp-2.5/dxi)**3. - &
               (4./3.)*pi*(radfp-3.5/dxi)**3. )/(1.*NL4)
       enddo
    enddo
    !
    return
  end subroutine coordsfp_interior
  !
end module mod_coordsfp
