! Runtime statistics in the layout of the CPL code (io.cpl, spanwShear
! variant), reduced on the device and written by the terminal rank.  With
! a Stokes layer, stokes_runtime.dat adds  time  energy_out  energy_in
! diss_out  diss_in  (region averages of <u_i u_i> and <grad u : grad u>,
! outside and inside |y - ly/2| < 8 delta, as the StokesLayer variant of
! io.cpl).
!
! Runtimedata columns:
!   time  meanflowx  meanflowy  S  S2  gamma_x  gamma_y  deltat  cfl*deltat
!   energy  diss  uw/2  vw/2
! variances_runtime.dat columns:
!   time  uu  vv  ww  uv
! damping_runtime.dat (hst_damping active):
!   time  theta_deg  ncol  energy_damped  power_damped
!   (ly/2 <u_i u_i> over the damped columns and ly/2 <2 g u_i u_i>, the
!   power the damping removes; zero for the hard mask)
! spectra/exz<i>.dat (with every field snapshot): the column energies
!   ix  iz  kx  kz  E  with E = ly/2 sum_y |u|^2 of the column (modes with
!   ix > 0 counted twice, as in energy), averaged over the dt_stat samples
!   since the previous snapshot.
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
  use hst_linsolve, only: line_solve, KIND_DY
  use hst_derivatives, only: s2_of, gamma_y_of
  use hst_stokes, only: stokes_active
  use hst_damping, only: damping_active, damping_theta, damping_soft_factor, ncol_damped

  implicit none
  private
  public :: outstats, open_runtimedata, close_runtimedata, spectrum_write

  integer, parameter :: unit_rt = 121, unit_var = 122, unit_sl = 123, unit_dp = 124
  real(C_DOUBLE), allocatable, save :: exz(:, :), exz_acc(:, :)   ! column energies: device scratch, host sum
  integer, save :: nacc = 0
  real(C_DOUBLE), save :: t_acc0 = 0.0d0

contains

  subroutine open_runtimedata()
    if (has_terminal) open (unit=unit_rt, file='Runtimedata', action='write', position='append')
    if (has_terminal) open (unit=unit_var, file='variances_runtime.dat', action='write', position='append')
    if (has_terminal .and. stokes_active()) open (unit=unit_sl, file='stokes_runtime.dat', action='write', position='append')
    if (has_terminal .and. damping_active()) open (unit=unit_dp, file='damping_runtime.dat', action='write', position='append')
    allocate (exz(-nz:nz, nx0:nxN), exz_acc(-nz:nz, nx0:nxN))
    exz = 0.0d0; exz_acc = 0.0d0; nacc = 0; t_acc0 = time
    !$omp target enter data map(to: exz)
  end subroutine open_runtimedata

  subroutine close_runtimedata()
    if (has_terminal) close (unit_rt)
    if (has_terminal) close (unit_var)
    if (has_terminal .and. stokes_active()) close (unit_sl)
    if (has_terminal .and. damping_active()) close (unit_dp)
    !$omp target exit data map(delete: exz)
    deallocate (exz, exz_acc)
  end subroutine close_runtimedata

  ! Writes the column energies averaged since the previous call (or the
  ! start) and restarts the average.
  subroutine spectrum_write(filename)
    character(len=*), intent(in) :: filename
    real(C_DOUBLE), allocatable :: full(:, :)
    integer(C_INT) :: ix, iz
    integer :: ierr, unit
    allocate (full(-nz:nz, 0:nx)); full = 0.0d0
    if (nacc > 0) full(:, nx0:nxN) = exz_acc/nacc
    call MPI_Allreduce(MPI_IN_PLACE, full, size(full), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    if (has_terminal) then
      open (newunit=unit, file=filename, action='write')
      write (unit, '(A,F12.5,A,F12.5,A,I0,A)') '# column energies ly/2 sum_y |u|^2 averaged over t = ', t_acc0, ' .. ', time, &
        ' (', nacc, ' samples);  ix  iz  kx  kz  E'
      do ix = 0, nx
        do iz = -nz, nz
          write (unit, '(I5,1X,I5,1X,ES15.7,1X,ES15.7,1X,ES15.7)') ix, iz, alfa0*ix, beta0*iz, 0.5d0*full(iz, ix)
        end do
      end do
      close (unit)
    end if
    deallocate (full)
    exz_acc = 0.0d0; nacc = 0; t_acc0 = time
  end subroutine spectrum_write

  subroutine outstats()
    integer(C_INT) :: ix, iy, iz, c
    real(C_DOUBLE) :: w, eps, uv, uu, vv, ww, vw, uw, grad, mfx, mfz
    real(C_DOUBLE) :: q_in, q_out, e_in, e_out, g_in, g_out, l_in, l_out
    real(C_DOUBLE) :: sums(16), q2, ed, pd, gsoft, qc
    complex(C_DOUBLE_COMPLEX) :: cu, cv, cw, dq
    integer :: ierr

    ! sums over modes and rows, each row weighted by its spacing dyl (so
    ! the results are integrals over the box height)
    eps = 0; uv = 0; uu = 0; vv = 0; ww = 0; vw = 0; uw = 0; mfx = 0; mfz = 0
    q_in = 0; q_out = 0; e_in = 0; e_out = 0; ed = 0; pd = 0
    gsoft = damping_soft_factor()
    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(V, k2, dyl, inlayer, gdamp, gsoft, nx0, nxN, nz, ny0, nyN) private(ix, iy, iz, w, cu, cv, cw, qc) &
    !$omp reduction(+:eps, uv, uu, vv, ww, vw, uw, mfx, mfz, q_in, q_out, e_in, e_out, ed, pd)
    do ix = nx0, nxN
      do iz = -nz, nz
        do iy = ny0, nyN
          cu = V(iy, iz, ix, 1); cv = V(iy, iz, ix, 2); cw = V(iy, iz, ix, 3)
          if (ix == 0 .and. iz == 0) then
            mfx = mfx + dreal(cu)*dyl(iy)
            mfz = mfz + dreal(cw)*dyl(iy)
            if (gdamp(iz, ix) > 0.0d0) then      ! damped mean profiles (damp_mean)
              qc = dyl(iy)*dreal(cu*conjg(cu) + cw*conjg(cw))
              ed = ed + qc; pd = pd + 2.0d0*gsoft*gdamp(iz, ix)*qc
            end if
            cycle
          end if
          w = 2.0d0*dyl(iy)
          if (ix == 0) w = dyl(iy)
          if (gdamp(iz, ix) > 0.0d0) then
            qc = w*dreal(cu*conjg(cu) + cv*conjg(cv) + cw*conjg(cw))
            ed = ed + qc; pd = pd + 2.0d0*gsoft*gdamp(iz, ix)*qc
          end if
          uu = uu + w*dreal(cu*conjg(cu))
          vv = vv + w*dreal(cv*conjg(cv))
          ww = ww + w*dreal(cw*conjg(cw))
          uv = uv + w*dreal(cu*conjg(cv))
          vw = vw + w*dreal(cw*conjg(cv))
          uw = uw + w*dreal(cu*conjg(cw))
          eps = eps + w*k2(iz, ix)*dreal(cu*conjg(cu) + cv*conjg(cv) + cw*conjg(cw))
          if (inlayer(iy) == 1) then
            q_in = q_in + w*dreal(cu*conjg(cu) + cv*conjg(cv) + cw*conjg(cw))
            e_in = e_in + w*k2(iz, ix)*dreal(cu*conjg(cu) + cv*conjg(cv) + cw*conjg(cw))
          else
            q_out = q_out + w*dreal(cu*conjg(cu) + cv*conjg(cv) + cw*conjg(cw))
            e_out = e_out + w*k2(iz, ix)*dreal(cu*conjg(cu) + cv*conjg(cv) + cw*conjg(cw))
          end if
        end do
      end do
    end do
    ! y derivatives of the three components, one at a time, into scratch
    do c = 1, 3
      call line_solve(KIND_DY, 0.0d0, V(:, :, :, c), rhs(:, :, :, 1))
      grad = 0; g_in = 0; g_out = 0
      !$omp target teams distribute parallel do collapse(3) default(none) &
      !$omp shared(rhs, dyl, inlayer, nx0, nxN, nz, ny0, nyN) private(ix, iy, iz, w, dq) reduction(+:grad, g_in, g_out)
      do ix = nx0, nxN
        do iz = -nz, nz
          do iy = ny0, nyN
            if (ix == 0 .and. iz == 0) cycle
            w = 2.0d0*dyl(iy)
            if (ix == 0) w = dyl(iy)
            dq = rhs(iy, iz, ix, 1)
            grad = grad + w*dreal(dq*conjg(dq))
            if (inlayer(iy) == 1) then
              g_in = g_in + w*dreal(dq*conjg(dq))
            else
              g_out = g_out + w*dreal(dq*conjg(dq))
            end if
          end do
        end do
      end do
      eps = eps + grad; e_in = e_in + g_in; e_out = e_out + g_out
    end do
    ! column energies (the 2-D spectrum), accumulated on the host
    !$omp target teams distribute parallel do collapse(2) default(none) &
    !$omp shared(V, dyl, exz, nx0, nxN, nz, ny0, nyN) private(ix, iy, iz, w, cu, cv, cw, qc)
    do ix = nx0, nxN
      do iz = -nz, nz
        qc = 0.0d0
        do iy = ny0, nyN
          cu = V(iy, iz, ix, 1); cv = V(iy, iz, ix, 2); cw = V(iy, iz, ix, 3)
          w = 2.0d0*dyl(iy)
          if (ix == 0) w = dyl(iy)
          qc = qc + w*dreal(cu*conjg(cu) + cv*conjg(cv) + cw*conjg(cw))
        end do
        exz(iz, ix) = qc
      end do
    end do
    !$omp target update from(exz)
    exz_acc = exz_acc + exz; nacc = nacc + 1
    ! over the ranks, then back into the named quantities
    sums = [eps, uv, uu, vv, ww, vw, uw, mfx, mfz, q_in, q_out, e_in, e_out, ed, pd, dble(ncol_damped)]
    call MPI_Allreduce(MPI_IN_PLACE, sums, 16, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    eps = sums(1); uv = sums(2); uu = sums(3); vv = sums(4); ww = sums(5); vw = sums(6); uw = sums(7)
    mfx = sums(8); mfz = sums(9); q_in = sums(10); q_out = sums(11); e_in = sums(12); e_out = sums(13)
    ed = sums(14); pd = sums(15)
    q2 = uu + vv + ww
    if (has_terminal) then
      ! energy = ly/2 <q2>, diss = ly/2 <grad u : grad u>, stresses ly/2 <..>, variances ly <..>
      write (*, '(F12.5,2X,ES11.4,2X,F8.4,4(2X,ES13.6))') time, deltat, cfl*deltat, &
        0.5d0*q2, 0.5d0*eps, 0.5d0*uv, 0.5d0*vw
      write (unit_rt, '(13(ES23.15,1X))') time, mfx, mfz, S, s2_of(time), S*time, gamma_y_of(time), deltat, cfl*deltat, &
        0.5d0*q2, 0.5d0*eps, 0.5d0*uv, 0.5d0*vw
      write (unit_var, '(5(ES23.15,1X))') time, uu, ww, vv, uw
      flush (unit_rt); flush (unit_var)
      if (stokes_active()) then
        ! region averages: q2 and grad u : grad u inside and outside |y - ly/2| < 8 delta (io.cpl)
        l_in = sum(dyl, mask=(inlayer == 1)); l_out = sum(dyl, mask=(inlayer == 0))
        write (unit_sl, '(5(ES23.15,1X))') time, q_out/l_out, q_in/l_in, e_out/l_out, e_in/l_in
        flush (unit_sl)
      end if
      if (damping_active()) then
        write (unit_dp, '(2(ES23.15,1X),I8,2(1X,ES23.15))') time, damping_theta(time)*180.0d0/3.141592653589793d0, &
          nint(sums(16)), 0.5d0*ed, 0.5d0*pd
        flush (unit_dp)
      end if
    end if
  end subroutine outstats

end module hst_stats
