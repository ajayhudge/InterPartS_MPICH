      module mod_param
      implicit none
!     Requirements for the number of processes in this parallellisation:
!     (i)   imax >= iu !for volume averaging
!     (ii)  jmax >= ju !for volume averaging
!     (iii) kmax must be an integer number of dims(1)*dims(2)
!     Other requirements:
!     (iv)  Largest prime number of itot and jtot should be 19 or lower

      integer,parameter :: ndims = 2
      integer,dimension(ndims),parameter :: dims = (/16,18/)

      integer itot,jtot,kmax
      !parameter(itot=432,jtot=432,kmax=864)
      parameter(itot=40,jtot=40,kmax=80)
      integer it1,jt1
      parameter(it1=itot+1,jt1=jtot+1)
      integer imax,jmax,i1,j1,k1
      parameter(imax=itot/dims(1),jmax=jtot/dims(2),i1=imax+1,j1=jmax+1,k1=kmax+1)

      real lx,ly,lz
      !parameter(lx=0.054,ly=0.054,lz=0.108)
      parameter(lx=0.0004,ly=0.0004,lz=0.0008)
      real dxi,dyi,dzi
      parameter(dxi=itot/lx,dyi=jtot/ly,dzi=kmax/lz)

      real radius,picon
      !parameter(radius=1e-3, picon=3.1415926535897932) 
      parameter(radius=5e-5, picon=3.1415926535897932)
!     Number of spheres
      integer NP
      parameter(NP=1) ! 1 x 1 x 1
      real retraction
      parameter(retraction=0.3/dxi) !lfps retracted from interface into interior of obstacle
!     Number of Lagrangian force point over surface sphere (at r=R-0.5/dxi !)
      integer NL
      parameter(NL =nint((picon/3.)*(12.*(((radius-retraction        )*dxi)**2) + 1. )))
!     nr lfp's of second shell
      integer NL2
      parameter(NL2=nint((picon/3.)*(12.*(((radius-retraction-1.2/dxi)*dxi)**2) + 1. )))
!     nr lfp's of third shell
      integer NL3
      parameter(NL3=nint((picon/3.)*(12.*(((radius-retraction-2.2/dxi)*dxi)**2) + 1. )))
!     nr lfp's of fourth shell
      integer NL4
      parameter(NL4=nint((picon/3.)*(12.*(((radius-retraction-3.2/dxi)*dxi)**2) + 1. )))
!     nr lfp's of 1st-4th shell
      integer NLtot
      parameter(NLtot=NL+NL2+NL3+NL4)

      real dt
      parameter(dt=0.01*radius**2/sqrt(1.*NL)) ! safety factor*R^2/sqrt(NL)

      real por
      parameter(por=1.-(1.*NP)*(4./3.)*picon*(radius**3)/(lx*ly*lz))

!     inversed Re-number
      real visc
      parameter(visc=1.62e-6)
!     inversed Re-number based on intrinsic bulk velocity and dimension cube

      real presgrad,fact1,fact2
      parameter(fact1 = (1.-(1.-por)**(1./3.))**3. )
      parameter(fact2 = (1.+(1.-por)**(1./3.)) )
      parameter(presgrad = visc*por*11.4*(1.-por)/(fact1*fact2) )  !for array of cubes!

      real gravacc
      parameter(gravacc=-9.81) !in units of u**2/l, working in the x-direction!
      real densratioconstant
      parameter(densratioconstant=1.0079) !rho_p/rho_f
      real momentofinertiaconstant
      parameter(momentofinertiaconstant=(2./5.)*(4./3.)*picon*(radius**3)*(radius**2)) !divided by mass density

!      real dt
!      parameter(dt=1.e-3)

      real rhof,usinf,stiffness,forcerange
      parameter(rhof=1.)  !fluid density normalized by fluid density
      parameter(usinf=1.) !terminal settling velocity normalized by bulk velocity
                          !Following Uhlmann (JCP, 2008), it is assumed that the terminal settling velocity
                          !is equal to the bulk velocity (=1). This can be approximately accomplished
                          !by controlling the value of the gravitational acceleration.
      parameter(stiffness=(8.e-4)*2.*radius/(rhof*(usinf**2)),forcerange=2./dxi)

      real uref,tref,lref
!      parameter(uref=sqrt(abs(gravacc*2.*radius)),tref=sqrt(abs(2.*radius/gravacc)),lref=lz)
      parameter(uref=1.,tref=lz,lref=lz)

!     Maximum number of particles on a process (for which a 
!     particular process is either master or slave)
      integer pmax
      parameter(pmax = NP)

      real offset
      parameter(offset = ( sqrt(3.*(1.5**2)) )/dxi + 0.01/dxi)
!     1.5/dxi is width of smoothed delta-function
!     Lagrangian force points located at r=radfp (assumed smaller or equal to radius) from centroid sphere
!     additional 0.01/dxi for safety


      end module mod_param
