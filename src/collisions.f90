module mod_collisions !VER O QUE ACONTECE SE TIVER NDT=40 E ITERMAX = 1
  use mpi
  use mod_common
  use mod_common_mpi
  implicit none
  private
  public collisions, lubrication
contains
  !
  subroutine collisions(p,q,qq,idp,idq,dist,deltax,deltay,deltaz, &
       u_nb,v_nb,w_nb,omx_nb,omy_nb,omz_nb,dtsub)
    implicit none
    integer, intent(in) :: p,q,qq,idp,idq
    real, intent(in) :: dist,deltax,deltay,deltaz, &
         u_nb,v_nb,w_nb,omx_nb,omy_nb,omz_nb, &
         dtsub
    integer :: i
    real :: distold,dtabs
    real :: nx,ny,nz,nxold,nyold,nzold,tx,ty,tz
    real :: dxn,dyn,dzn,dxt,dyt,dzt
    real :: dun,dvn,dwn,dut,dvt,dwt,vabs,vtabs
    real :: fxn,fyn,fzn,fnabs,fxt,fyt,fzt,ftabs
    real :: psi,psi_crit
    real :: deltan,kn,etan,kt,etat,muc
    real :: torqx,torqy,torqz
    real :: hx,hy,hz,aa,bb,cc
    real :: h11,h12,h13, &
         h21,h22,h23, &
         h31,h32,h33
    real :: meffn,mefft
    !
    ! vector (nx,ny,nz) is normal unit vector, pointing from particle p towards particle q
    !
    nx = deltax/dist
    ny = deltay/dist
    nz = deltaz/dist
    !
    if (idq.le.np) then
       deltan = 2*radius-dist
       meffn = meffn_ss
       mefft = mefft_ss
       kn = kn_ss
       kt = kt_ss
       etan = etan_ss
       etat = etat_ss
       muc = muc_ss
       psi_crit = psi_crit_ss
    else
       deltan = radius-dist
       meffn = meffn_sw
       mefft = mefft_sw
       kn = kn_sw
       kt = kt_sw
       etan = etan_sw
       etat = etat_sw
       muc = muc_sw
       psi_crit = psi_crit_sw
    endif
    !
    ! relative velocity between particle p and q (Nb: p - q !)
    !
    dun = (ap(p)%u+(ap(p)%omy*radius*nz) - (ap(p)%omz*radius*ny)) - &
         (u_nb + (omy_nb*radius*(-nz)) - (omz_nb*radius*(-ny)))
    dvn = (ap(p)%v-(ap(p)%omx*radius*nz) + (ap(p)%omz*radius*nx)) - &
         (v_nb - (omx_nb*radius*(-nz)) + (omz_nb*radius*(-nx)))
    dwn = (ap(p)%w+(ap(p)%omx*radius*ny) - (ap(p)%omy*radius*nx)) - &
         (w_nb + (omx_nb*radius*(-ny)) - (omy_nb*radius*(-nx)))
    vabs = dun*nx + dvn*ny + dwn*nz
    !
    ! computation of contact forces based on soft-sphere model
    !
    fxn = - kn*(deltan*nx) - etan*(vabs*nx)
    fyn = - kn*(deltan*ny) - etan*(vabs*ny)
    fzn = - kn*(deltan*nz) - etan*(vabs*nz)
    !
    ! relative tangential velocity
    !
    dut = dun - vabs*nx
    dvt = dvn - vabs*ny
    dwt = dwn - vabs*nz
    !
    vtabs = sqrt(dut**2.+dvt**2.+dwt**2.)
    !
    ! prediction of contact instant skipped for now for simplicity's sake
    !
    distold = sqrt(op(p)%dx(qq)**2.+op(p)%dy(qq)**2.+op(p)%dz(qq)**2.)
    if(distold.eq.0) then
       nxold = nx
       nyold = ny
       nzold = nz
    else
       nxold = op(p)%dx(qq)/distold
       nyold = op(p)%dy(qq)/distold
       nzold = op(p)%dz(qq)/distold
    endif
    !
    ! computation of rotation matrix hij
    !
    if(ap(p)%firstc(qq).eq.idq) then
       psi = ap(p)%psi(qq)
       hx = ny*nzold-nyold*nz
       hy = nz*nxold-nzold*nx
       hz = nx*nyold-nxold*ny
       aa = sqrt(hx**2. + hy**2. + hz **2.)
       cc = cos(asin(aa))
       bb = 1. - cc
       if(aa.eq.0.) then
          hx = 0.
          hy = 0.
          hz = 0. 
       else
          hx = hx/aa
          hy = hy/aa
          hz = hz/aa
       endif
       h11 = bb*hx**2. + cc
       h12 = bb*hx*hy - aa*hz
       h13 = bb*hx*hz + aa*hy
       h21 = bb*hx*hy + aa*hz
       h22 = bb*hy**2. + cc
       h23 = bb*hy*hz - aa*hx
       h31 = bb*hx*hz - aa*hy
       h32 = bb*hy*hz + aa*hx
       h33 = bb*hz**2. + cc
       !  
       ! tangential displacement evolved in time
       !
       dxt = op(p)%dxt(qq)*h11 + op(p)%dyt(qq)*h21 + op(p)%dzt(qq)*h31 + &
            0.5*dtsub*rkcoeffab(rkiter)*(dut + &
            op(p)%dut(qq)*h11 + op(p)%dvt(qq)*h21 + op(p)%dwt(qq)*h31)
       dyt = op(p)%dxt(qq)*h12 + op(p)%dyt(qq)*h22 + op(p)%dzt(qq)*h32 + &
            0.5*dtsub*rkcoeffab(rkiter)*(dvt + &
            op(p)%dut(qq)*h12 + op(p)%dvt(qq)*h22 + op(p)%dwt(qq)*h32)
       dzt = op(p)%dxt(qq)*h13 + op(p)%dyt(qq)*h23 + op(p)%dzt(qq)*h33 + &
            0.5*dtsub*rkcoeffab(rkiter)*(dwt + &
            op(p)%dut(qq)*h13 + op(p)%dvt(qq)*h23 + op(p)%dwt(qq)*h33)
    else
       !  ap(p)%firstc(qq) = idq
       psi = vtabs/abs(vabs)
       ap(p)%psi(qq) = psi
       dxt = 0.
       dyt = 0.
       dzt = 0.
    endif
    dtabs = sqrt(dxt**2.+dyt**2.+dzt**2.)
    !
    ! computation of tangential force
    !
    !if(psi.lt.psi_crit) then ! stick
    !  fxt = - (kt*dxt) - (etat*dut)
    !  fyt = - (kt*dyt) - (etat*dvt)
    !  fzt = - (kt*dzt) - (etat*dwt)
    !  print*,'stick'
    !else ! slip
    !  fnabs = sqrt(fxn**2.+fyn**2.+fzn**2.)
    !  fxt = - muc*fnabs*tx
    !  fyt = - muc*fnabs*ty
    !  fzt = - muc*fnabs*tz
    !  print*,'slip'
    !endif
    !
    fxt = - (kt*dxt) - (etat*dut)
    fyt = - (kt*dyt) - (etat*dvt)
    fzt = - (kt*dzt) - (etat*dwt)
    fnabs = sqrt(fxn**2.+fyn**2.+fzn**2.)
    ftabs = sqrt(fxt**2.+fyt**2.+fzt**2.)
    tx = - fxt/ftabs
    ty = - fyt/ftabs
    tz = - fzt/ftabs
    if(ftabs.le.muc*fnabs) then
       !  write(*,*) 'Sticking'
    else
       dxt = muc*fnabs*tx/kt
       dyt = muc*fnabs*ty/kt
       dzt = muc*fnabs*tz/kt
       fxt = - (kt*dxt) - (etat*dut)
       fyt = - (kt*dyt) - (etat*dvt)
       fzt = - (kt*dzt) - (etat*dwt)
       ftabs = sqrt(fxt**2.+fyt**2.+fzt**2.)
       if(ftabs.le.muc*fnabs) then
          !we were lucky
          !    write(*,*) 'frontier between stick and slip'
       else
          fxt = -muc*fnabs*tx
          fyt = -muc*fnabs*ty
          fzt = -muc*fnabs*tz
          !    write(*,*) 'Slipping'
       endif
    endif
    !
    ! computation of collision torque
    !
    torqx = radius*(ny*fzt-nz*fyt)
    torqy = radius*(nz*fxt-nx*fzt)
    torqz = radius*(nx*fyt-ny*fxt)
    !
    ! add contribution from this contact to the total contact forces/torques
    !
    ap(p)%colfx = ap(p)%colfx + fxn! + fxt*0.
    ap(p)%colfy = ap(p)%colfy + fyn! + fyt*0.
    ap(p)%colfz = ap(p)%colfz + fzn! + fzt*0.
    !
    ap(p)%coltx = ap(p)%coltx! + torqx*0.
    ap(p)%colty = ap(p)%colty! + torqy*0.
    ap(p)%coltz = ap(p)%coltz! + torqz*0.
    !
    ! update variables needed for integrating the tangential displacement
    !
    ap(p)%dx(qq) = deltax
    ap(p)%dy(qq) = deltay
    ap(p)%dz(qq) = deltaz
    ap(p)%dxt(qq) = dxt
    ap(p)%dyt(qq) = dyt
    ap(p)%dzt(qq) = dzt
    ap(p)%dut(qq) = dut
    ap(p)%dvt(qq) = dvt
    ap(p)%dwt(qq) = dwt
    !
    return
  end subroutine collisions
  !
  subroutine lubrication(p,q,idp,idq,nx,ny,nz,eps,u_nb,v_nb,w_nb,omx_nb,omy_nb,omz_nb)
    implicit none
    integer, intent(in) :: p,q,idp,idq
    real, intent(in) :: nx,ny,nz,eps, &
         u_nb,v_nb,w_nb, &
         omx_nb,omy_nb,omz_nb
    !
    ! 1-> squeezing direction; 2 and 3 -> shearing directions
    ! see paper by dance and Maxey
    !
    real :: a11,a22,a33,b23,b32,c23,c32,d11,d22,d33           
    real :: a11a,a22a,a33a,b23a,b32a,c23a,c32a,d11a,d22a,d33a 
    real :: a11b,a22b,a33b,b23b,b32b,c23b,c32b,d11b,d22b,d33b 
    real :: u1a,u2a,u3a,omg1a,omg2a,omg3a 
    real :: u1b,u2b,u3b,omg1b,omg2b,omg3b
    real :: t1x,t1y,t1z,t2x,t2y,t2z
    real :: f1,f2,f3,t1,t2,t3
    real :: dun,dvn,dwn,dut,dvt,dwt,vabs,vtabs
    !
    ! determine the unit vectors t1x and t2x from nx,ny,nz
    ! and the velocity differences
    !
    dun = (ap(p)%u+(ap(p)%omy*radius*nz) - (ap(p)%omz*radius*ny)) - &
         (u_nb + (omy_nb*radius*(-nz)) - (omz_nb*radius*(-ny)))
    dvn = (ap(p)%v-(ap(p)%omx*radius*nz) + (ap(p)%omz*radius*nx)) - &
         (v_nb - (omx_nb*radius*(-nz)) + (omz_nb*radius*(-nx)))
    dwn = (ap(p)%w+(ap(p)%omx*radius*ny) - (ap(p)%omy*radius*nx)) - &
         (w_nb + (omx_nb*radius*(-ny)) - (omy_nb*radius*(-nx)))
    vabs = dun*nx + dvn*ny + dwn*nz
    dut = dun - vabs*nx
    dvt = dvn - vabs*ny
    dwt = dwn - vabs*nz
    !
    vtabs = sqrt(dut**2.+dvt**2.+dwt**2.)
    !
    ! t1i has the direction of the tangential velocity at contact point
    !
    t1x = dut/vtabs
    t1y = dvt/vtabs
    t1z = dwt/vtabs
    !
    ! t2i results from the cross product of ni and t1i
    !
    t2x = ny*t1z-nz*t1y
    t2y = nz*t1x-nx*t1z
    t2z = nx*t1y-ny*t1x
    !
    ! compute velocity differences in local reference frame
    !
    u1a = ap(p)%u*nx + ap(p)%v*ny + ap(p)%w*nz
    u1b = u_nb*nx + v_nb*ny + w_nb*nz
    u2a = ap(p)%u*t1x + ap(p)%v*t1y + ap(p)%w*t1z
    u2b = u_nb*t1x + v_nb*t1y + w_nb*t1z
    u3a = ap(p)%u*t2x + ap(p)%v*t2y + ap(p)%w*t2z
    u3b = u_nb*t2x + v_nb*t2y + w_nb*t2z
    omg1a = ap(p)%omx*nx + ap(p)%omy*ny+ap(p)%omz*nz
    omg1b = omx_nb*nx + omy_nb*ny+omz_nb*nz
    omg2a = ap(p)%omx*t1x + ap(p)%omy*t1y+ap(p)%omz*t1z
    omg2b = omx_nb*t1x + omy_nb*t1y+omz_nb*t1z
    omg3a = ap(p)%omx*t2x + ap(p)%omy*t2y+ap(p)%omz*t2z
    omg3b = omx_nb*t2x + omy_nb*t2y+omz_nb*t2z
    !
    ! compute lubrication forces
    !
    if(idq.gt.np) then ! particle-wall collision
       if(eps.lt.eps_sat_pw) then
          a11 = a11_sat_pw
          a22 = a22_sat_pw
          a33 = a33_sat_pw
          b23 = b23_sat_pw
          b32 = b32_sat_pw
          c23 = c23_sat_pw
          c32 = c32_sat_pw
          d11 = d11_sat_pw
          d22 = d22_sat_pw
          d33 = d33_sat_pw
       else
          a11 = -1./eps+1./5.*log(eps)+1./21.*eps*log(eps)-.9713
          a22 = 8./15.*log(eps)+64./375.*eps*log(eps)-.952
          a33 = a22
          b23 = -2./15.*log(eps)-86./375.*eps*log(eps)-.257
          b32 = -b23
          c23 = b32
          c32 = b23
          d11 = 1./2.*eps*log(eps)-1.202
          d22 = 2./5.*log(eps)+66./125.*eps*log(eps)-.371
          d33 = 2./5.*log(eps)+66./125.*eps*log(eps)-.371
       endif
       f1 = (a11*u1a) ! squeezing
       f2 = (a22*u2a + radius*b23*omg3a) ! translational + rotational shearing
       f3 = (a33*u3a + radius*b32*omg2a) ! translational + rotational shearing
       t1 = (radius*d11*omg1a) ! torsoidal shearing
       t2 = (c23*u3a+radius*d22*omg2a) ! translational + rotational shearing
       t3 = (c32*u2a+radius*d33*omg3a) ! translational + rotational shearing
    else
       if(eps.lt.eps_sat_pp) then
          a11  = a11_sat_pp
          a11a = a11a_sat_pp
          a11b = a11b_sat_pp
          a22  = a22_sat_pp
          a22a = a22a_sat_pp
          a22b = a22b_sat_pp
          a33a = a33a_sat_pp
          a33b = a33b_sat_pp
          b23  = b23_sat_pp
          b23a = b23a_sat_pp
          b23b = b23b_sat_pp
          b32a = b32a_sat_pp
          b32b = b32b_sat_pp
          c23a = c23a_sat_pp
          c23b = c23b_sat_pp
          c32a = c32a_sat_pp
          c32b = c32b_sat_pp
          d11  = d11_sat_pp
          d11a = d11a_sat_pp
          d11a = d11b_sat_pp
          d22a = d22a_sat_pp
          d22b = d22b_sat_pp
          d33a = d33a_sat_pp
          d33b = d33b_sat_pp
       else
          a11  = -1./4.*eps**(-1.)+9./40.*log(eps)+3./112.*eps*log(eps)
          a11a = a11-.995
          a11b = -a11+.350
          a22  = 1./6.*log(eps)
          a22a = a22-.998
          a22b = -a22+.274
          a33a = a22a
          a33b = a22b
          b23  = -1./6.*log(eps)-1./12.*eps*log(eps)
          b23a = b23-.159
          b23b = -b23+.001
          b32a = -b23a
          b32b = -b23b
          c23a = b32a
          c23b = b32b
          c32a = b23a
          c32b = b23b
          d11  = 1./8.*eps*log(eps)
          d11a = d11
          d11b = -d11
          d22a = 1./5.*log(eps)+47./250.*eps*log(eps)-.703
          d22b = -1./20.*log(eps)+31./500.*eps*log(eps)-.027
          d33a = 1./5.*log(eps)+47./250.*eps*log(eps)-.703
          d33b = -1./20.*log(eps)+31./500.*eps*log(eps)-.027
       endif
       f1 = (a11a*u1a+a11b*u1b) ! squeezing 
       f2 = (a22a*u2a+a22b*u2b+radius*(b23a*omg3a-b23b*omg3b)) ! translational + rotational shearing
       f3 = (a33a*u3a+a33b*u3b+radius*(b32a*omg2a-b32b*omg2b)) ! translational + rotational shearing
       t1 = (radius*(d11a*omg1a-d11b*omg1b)) ! torsoidal shearing
       t2 = (c23a*u3a+c23b*u3b+radius*(d22a*omg2a-d22b*omg2b)) ! translational + rotational shearing
       t3 = (c32a*u2a+c32b*u2b+radius*(d33a*omg3a-d33b*omg3b)) ! translational + rotational shearing
    endif
    !
    if(idq.gt.np) then ! particle-wall collision
       a11 = a11_ini_pw
       a22 = a22_ini_pw
       a33 = a33_ini_pw
       b23 = b23_ini_pw
       b32 = b32_ini_pw
       c23 = c23_ini_pw
       c32 = c32_ini_pw
       d11 = d11_ini_pw
       d22 = d22_ini_pw
       d33 = d33_ini_pw
       f1 = coeff_f*(f1-(a11*u1a))
       f2 = coeff_f*(f2-(a22*u2a + radius*b23*omg3a))
       f3 = coeff_f*(f3-(a33*u3a + radius*b32*omg2a))
       t1 = coeff_t*(t1-(radius*d11*omg1a))
       t2 = coeff_t*(t2-(c23*u3a+radius*d22*omg2a))
       t3 = coeff_t*(t3-(c32*u2a+radius*d33*omg3a))
    else
       a11  = a11_ini_pp
       a11a = a11a_ini_pp
       a11b = a11b_ini_pp
       a22  = a22_ini_pp
       a22a = a22a_ini_pp
       a22b = a22b_ini_pp
       a33a = a33a_ini_pp
       a33b = a33b_ini_pp
       b23  = b23_ini_pp
       b23a = b23a_ini_pp
       b23b = b23b_ini_pp
       b32a = b32a_ini_pp
       b32b = b32b_ini_pp
       c23a = c23a_ini_pp
       c23b = c23b_ini_pp
       c32a = c32a_ini_pp
       c32b = c32b_ini_pp
       d11  = d11_ini_pp
       d11a  = d11a_ini_pp
       d11b  = d11b_ini_pp
       d22a = d22a_ini_pp
       d22b = d22b_ini_pp
       d33a = d33a_ini_pp
       d33b = d33b_ini_pp
       f1 = coeff_f*(f1-(a11a*u1a+a11b*u1b))
       f2 = coeff_f*(f2-(a22a*u2a+a22b*u2b+radius*(b23a*omg3a-b23b*omg3b)))
       f3 = coeff_f*(f3-(a33a*u3a+a33b*u3b+radius*(b32a*omg2a-b32b*omg2b)))
       t1 = coeff_t*(t1-(radius*(d11a*omg1a-d11b*omg1b)))
       t2 = coeff_t*(t2-(c23a*u3a+c23b*u3b+radius*(d22a*omg2a-d22b*omg2b)))
       t3 = coeff_t*(t3-(c32a*u2a+c32b*u2b+radius*(d33a*omg3a-d33b*omg3b)))
       !
    endif
    !
    ! project forces in global reference frame
    ! non-normal interactions neglected for now
    !
    ap(p)%colfx = ap(p)%colfx + f1*nx! + f2*t1x + f3*t2x
    ap(p)%colfy = ap(p)%colfy + f1*ny! + f2*t1y + f3*t2y
    ap(p)%colfz = ap(p)%colfz + f1*nz! + f2*t1z + f3*t2z
    !ap(p)%coltx = ap(p)%coltx + t1*nx + t2*t1x + t3*t2x
    !ap(p)%colty = ap(p)%coltx + t1*ny + t2*t1y + t3*t2y
    !ap(p)%coltz = ap(p)%coltz + t1*nz + t2*t1z + t3*t2z
    !
    return
  end subroutine lubrication
  !
end module mod_collisions
