module mod_chkdiv
  implicit none
  private
  public chkdiv
contains
  !
  subroutine chkdiv
    use mpi
    use mod_param
    use mod_common
    use mod_common_mpi
    implicit none
    real :: div,divtot,divtot_all,divmax(2),divmax_all(2)
    integer :: i,j,k,im,jm,km
    !
    divmax(1) = 0.
    divmax(2) = 1.*myid
    divtot = 0.
    im = 0
    jm = 0
    km = 0
    !
    do k=1,kmax
       do j=1,jmax
          do i=1,imax
             div = (wnew(i,j,k)-wnew(i,j,k-1))*dzi + &
                   (vnew(i,j,k)-vnew(i,j-1,k))*dyi + &
                   (unew(i,j,k)-unew(i-1,j,k))*dxi
             divtot = divtot+div
             div = abs(div)
             if(div.gt.divmax(1)) then
                divmax(1) = div ! update local max divergence (SS)
                im = i
                jm = j
                km = k
             endif
          enddo
       enddo
    enddo
    !
    call mpi_allreduce(divtot,divtot_all,1,mpi_real8,mpi_sum,MPI_COMM_WORLD,error) ! global total divergence (SS)
    call mpi_allreduce(divmax,divmax_all,1,mpi_2double_precision,mpi_maxloc,MPI_COMM_WORLD,error) ! global max divergence and rank (SS)
    !
    if (myid.eq.int(divmax_all(2))) then
       write(6,111) zstart(1)-1+im, zstart(2)-1+jm,km
       write(6,222) divtot_all,divmax_all(1),int(divmax_all(2))
111    format('Maximal divergence at i = ',I5,' j = ', I5,' k = ',I5)
222    format('Divergence: Tot = ',e13.6,' Max = ',e13.6,' Rank = ',I3)
    endif
    !
    if (divtot_all .gt. 1.e-3) then
       if (myid .eq. divmax_all(2)) then
          write(6,*) 'Fatal error: total divergence > 1.e-3!'
          write(6,*) 'Program aborted...'
       endif
       call mpi_finalize(error)
       stop
    endif
    !
    return
  end subroutine chkdiv
  !
end module mod_chkdiv
