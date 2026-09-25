! Runtime statistics in the layout of the CPL code (io.cpl, spanwShear
! variant), reduced on the device and written by the terminal rank.
!
! Runtimedata columns:
!   time  meanflowx  meanflowy  S  S2  gamma_x  gamma_y  deltat  cfl*deltat
!   energy  diss  uw/2  vw/2
! variances_runtime.dat columns:
!   time  uu  vv  ww  uv
! all in CPL naming (their v is our spanwise w, their w our vertical v) and
! as integrals over the box height ly: energy = ly/2 <u_i u_i>, diss =
! ly/2 <du_i/dx_j du_i/dx_j> (not multiplied by nu; here from the compact
! derivatives, in CPL from centred differences), uw/2 = ly/2 <u v>,
! vw/2 = ly/2 <w v>, and the variances uu = ly <u u> etc.  meanflowx/y are
! the integrals of the mean profiles.  The (0,0) mode is excluded from the
! fluctuation statistics; modes with ix > 0 count twice.
module hst_stats

  use, intrinsic :: iso_c_binding
  use mpi_f08
  use hst_params
  use hst_linsolve, only: apply_dy

  implicit none
  private
  public :: outstats, open_runtimedata, close_runtimedata

  integer, parameter :: unit_rt = 121, unit_var = 122

contains

  subroutine open_runtimedata()
    if (has_terminal) open (unit=unit_rt, file='Runtimedata', action='write', position='append')
    if (has_terminal) open (unit=unit_var, file='variances_runtime.dat', action='write', position='append')
  end subroutine open_runtimedata

  subroutine close_runtimedata()
    if (has_terminal) close (unit_rt)
    if (has_terminal) close (unit_var)
  end subroutine close_runtimedata

  subroutine outstats()
    integer(C_INT) :: ix, iy, iz, c
    real(C_DOUBLE) :: w, eps, uv, uu, vv, ww, vw, uw, grad, mfx, mfz
    real(C_DOUBLE) :: sums(9), glob(9)      ! (not s/S: Fortran is case-insensitive and S is the shear)
    complex(C_DOUBLE_COMPLEX) :: cu, cv, cw, dq
    integer :: ierr

    eps = 0; uv = 0; uu = 0; vv = 0; ww = 0; vw = 0; uw = 0; mfx = 0; mfz = 0
    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(V, k2, nx0, nxN, nz, ny) private(ix, iy, iz, w, cu, cv, cw) &
    !$omp reduction(+:eps, uv, uu, vv, ww, vw, uw, mfx, mfz)
    do ix = nx0, nxN
      do iz = -nz, nz
        do iy = 0, ny - 1
          cu = V(iy, iz, ix, 1); cv = V(iy, iz, ix, 2); cw = V(iy, iz, ix, 3)
          if (ix == 0 .and. iz == 0) then
            mfx = mfx + dreal(cu)
            mfz = mfz + dreal(cw)
            cycle
          end if
          w = 2.0d0
          if (ix == 0) w = 1.0d0
          uu = uu + w*dreal(cu*conjg(cu))
          vv = vv + w*dreal(cv*conjg(cv))
          ww = ww + w*dreal(cw*conjg(cw))
          uv = uv + w*dreal(cu*conjg(cv))
          vw = vw + w*dreal(cw*conjg(cv))
          uw = uw + w*dreal(cu*conjg(cw))
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
    ! integrals over the box height: sum * dy = sum * ly/ny
    sums = [uu + vv + ww, eps, uv, vw, uu, vv, ww, uw, mfz]*(ly/ny)
    glob = 0
    call MPI_Allreduce(sums, glob, 9, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    mfx = mfx*(ly/ny)
    call MPI_Allreduce(MPI_IN_PLACE, mfx, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    if (has_terminal) then
      ! energy = ly/2 <q2>, diss = ly/2 <grad u : grad u>, stresses ly/2 <..>, variances ly <..>
      write (*, '(F12.5,2X,ES11.4,2X,F8.4,4(2X,ES13.6))') time, deltat, cfl*deltat, &
        0.5d0*glob(1), 0.5d0*glob(2), 0.5d0*glob(3), 0.5d0*glob(4)
      write (unit_rt, '(13(ES23.15,1X))') time, mfx, glob(9), S, 0.0d0, S*time, 0.0d0, deltat, cfl*deltat, &
        0.5d0*glob(1), 0.5d0*glob(2), 0.5d0*glob(3), 0.5d0*glob(4)
      write (unit_var, '(5(ES23.15,1X))') time, glob(5), glob(7), glob(6), glob(8)
      flush (unit_rt); flush (unit_var)
    end if
  end subroutine outstats

end module hst_stats
