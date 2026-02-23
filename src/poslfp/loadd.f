      module mod_loadd
      use mod_param
      implicit none
      contains
!---------------------------------------------------------------------

      subroutine loadd(in,nr,phi,theta)
!
!  BINARY I/O ROUTINE
!-------------------------------------------------------------
!-------------------------------------------------------------
!
      implicit none
      integer l
      integer in,reclengte,nr
!      parameter(reclengte = 8)
      real phi(:),theta(:)
      ! wwvv
      real xdum
      inquire (iolength=reclengte) xdum
      ! /wwvv

      if (in.eq.0) then
        open(15,file='stopiter',access='direct',
     &       recl=(2*NL+1)*reclengte)
        read(15,rec=1)
     &  (theta(l),l=1,NL),
     &  (phi(l),l=1,NL),
     &  nr
        close(15)
      endif
      if (in.eq.1) then
        open(25,file='stopiter',access='direct',
     &       recl=(2*NL+1)*reclengte)
        write(25,rec=1)
     &  (theta(l),l=1,NL),
     &  (phi(l),l=1,NL),
     &  nr
        close(25)
      endif

      return
      end subroutine loadd
      end module mod_loadd
