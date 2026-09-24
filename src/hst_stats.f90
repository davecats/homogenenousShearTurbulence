! Runtime statistics: box averages of the fluctuation field, reduced on the
! device and written to Runtimedata by the terminal rank.
!
! Columns of Runtimedata:
!   time  deltat  cfl  q2  eps  uv  uu  vv  ww
! with q2 = <u_i u_i> (twice the kinetic energy), eps = ni <du_i/dx_j du_i/dx_j>
! (y derivatives with the compact scheme, through apply_dy), and the
! Reynolds stresses.  The (0,0) mode -- the mean profiles -- is excluded.
! Modes with ix > 0 count twice (their conjugates are not stored).
module hst_stats

  use, intrinsic :: iso_c_binding
  use mpi_f08
  use hst_params
  use hst_linsolve, only: apply_dy

  implicit none
  private
  public :: outstats, open_runtimedata, close_runtimedata

  integer, parameter :: unit_rt = 121

contains

  subroutine open_runtimedata()
    if (has_terminal) open (unit=unit_rt, file='Runtimedata', action='write', position='append')
  end subroutine open_runtimedata

  subroutine close_runtimedata()
    if (has_terminal) close (unit_rt)
  end subroutine close_runtimedata

  subroutine outstats()
    integer(C_INT) :: ix, iy, iz, c
    real(C_DOUBLE) :: w, eps, uv, uu, vv, ww, grad
    real(C_DOUBLE) :: s(6), g(6)
    complex(C_DOUBLE_COMPLEX) :: cu, cv, cw, dq
    integer :: ierr

    eps = 0; uv = 0; uu = 0; vv = 0; ww = 0
    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(V, k2, nx0, nxN, nz, ny) private(ix, iy, iz, w, cu, cv, cw) reduction(+:eps, uv, uu, vv, ww)
    do ix = nx0, nxN
      do iz = -nz, nz
        do iy = 0, ny - 1
          if (ix == 0 .and. iz == 0) cycle
          w = 2.0d0
          if (ix == 0) w = 1.0d0
          cu = V(iy, iz, ix, 1); cv = V(iy, iz, ix, 2); cw = V(iy, iz, ix, 3)
          uu = uu + w*dreal(cu*conjg(cu))
          vv = vv + w*dreal(cv*conjg(cv))
          ww = ww + w*dreal(cw*conjg(cw))
          uv = uv + w*dreal(cu*conjg(cv))
          eps = eps + w*k2(iz, ix)*dreal(cu*conjg(cu) + cv*conjg(cv) + cw*conjg(cw))
        end do
      end do
    end do
    ! y derivatives of the three components, one at a time, into scratch
    do c = 1, 3
      call apply_dy(c, memrhs(:, :, :, 1))
      grad = 0
      !$omp target teams distribute parallel do collapse(3) default(none) &
      !$omp shared(memrhs, nx0, nxN, nz, ny) private(ix, iy, iz, w, dq) reduction(+:grad)
      do ix = nx0, nxN
        do iz = -nz, nz
          do iy = 0, ny - 1
            if (ix == 0 .and. iz == 0) cycle
            w = 2.0d0
            if (ix == 0) w = 1.0d0
            dq = memrhs(iy, iz, ix, 1)
            grad = grad + w*dreal(dq*conjg(dq))
          end do
        end do
      end do
      eps = eps + grad
    end do
    eps = ni*eps
    s = [uu + vv + ww, eps, uv, uu, vv, ww]/ny
    call MPI_Allreduce(s, g, 6, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    if (has_terminal) then
      write (*, '(F12.5,2X,ES11.4,2X,F8.4,6(2X,ES13.6))') time, deltat, cfl*deltat, g
      write (unit_rt, '(F14.7,2X,ES14.7,2X,F10.6,6(2X,ES16.9))') time, deltat, cfl*deltat, g
      flush (unit_rt)
    end if
  end subroutine outstats

end module hst_stats
