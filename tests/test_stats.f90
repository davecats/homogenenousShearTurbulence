! Averaging operator of the runtime statistics (hst_stats).
!
!   mpirun -np N build-cpu/test_stats [hst.in]
!
! The fluctuation is the field minus its box mean, the (0,0,0) mode only: the
! y-dependent part of the (0,0) column (the plane-averaged profiles u00(y),
! w00(y)) belongs to it, as in the post-processor (post/cpl, mod_rey_decomp.cpl).
! The random start has no (0,0) mode; this test puts profiles with a large box
! mean into it, runs outstats and spectrum_write once, and compares with a
! reference built the other way round: the box mean is subtracted from the
! field itself and then every mode is summed with no exception, the y
! derivatives from a line solve of the subtracted field.  Compared: uu, vv, ww,
! uw (variances_runtime.dat), energy, diss, uv, vw (Runtimedata), and the (0,0)
! column energy of the spectrum file.  The same numbers with the old operator
! (the whole (0,0) column dropped) must be different, or the test would not
! see the change.
program test_stats

  use, intrinsic :: iso_c_binding
  use mpi_f08
  use hst_params
  use test_common
  use hst_mpi
  use hst_derivatives
  use hst_linsolve
  use hst_initial
  use hst_io, only: make_output_dirs
  use hst_stats

  implicit none

  real(C_DOUBLE), parameter :: PI = 3.141592653589793d0
  integer :: ierr, m, ix, iy, iz, c, unit, ios
  integer(C_INT) :: ixx, izz
  real(C_DOUBLE) :: loc(2), ubar, wbar, lyw, w, grad, kx, kz
  real(C_DOUBLE) :: sum_new(8), sum_old(8), e00, got_rt(13), got_var(5), got_e00, ref(9), refold(9), errmax, d_old
  complex(C_DOUBLE_COMPLEX) :: cu, cv, cw, dq
  character(len=256) :: line
  logical :: passed

  call test_start()
  call test_setup()
  time = 0.0d0
  call generate_initial_field()

  ! profiles into the (0,0) column: a large box mean and y-dependent parts
  if (has_average) then
    do iy = ny0 - 2, nyN + 2
      V(iy, 0, 0, 1) = dcmplx(0.30d0 + 0.20d0*sin(2*PI*y(iy)/ly) + 0.10d0*cos(4*PI*y(iy)/ly), 0.0d0)
      V(iy, 0, 0, 3) = dcmplx(-0.15d0 + 0.25d0*cos(2*PI*y(iy)/ly) + 0.05d0*sin(6*PI*y(iy)/ly), 0.0d0)
      V(iy, 0, 0, 2) = 0.0d0
    end do
  end if
  !$omp target update to(V)
  do m = 1, 3
    call fill_ghosts(m)
  end do

  ! the code under test
  call make_output_dirs()
  call open_runtimedata()
  call outstats()
  call spectrum_write('exz_test.dat')
  call close_runtimedata()

  ! reference: subtract the box mean of u and w from the field itself
  lyw = sum(dyl)
  loc = 0.0d0
  if (has_average) then
    do iy = ny0, nyN
      loc(1) = loc(1) + dreal(V(iy, 0, 0, 1))*dyl(iy)
      loc(2) = loc(2) + dreal(V(iy, 0, 0, 3))*dyl(iy)
    end do
  end if
  call MPI_Allreduce(MPI_IN_PLACE, loc, 2, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  ubar = loc(1)/lyw; wbar = loc(2)/lyw
  if (has_average) then
    V(:, 0, 0, 1) = V(:, 0, 0, 1) - ubar
    V(:, 0, 0, 3) = V(:, 0, 0, 3) - wbar
  end if
  !$omp target update to(V)
  do m = 1, 3
    call fill_ghosts(m)
  end do

  ! order: uu vv ww uv vw uw eps(spectral+y) q; and the (0,0) column energy
  sum_new = 0.0d0; sum_old = 0.0d0; e00 = 0.0d0
  do ix = nx0, nxN
    do iz = -nz, nz
      do iy = ny0, nyN
        cu = V(iy, iz, ix, 1); cv = V(iy, iz, ix, 2); cw = V(iy, iz, ix, 3)
        w = 2.0d0*dyl(iy); if (ix == 0) w = dyl(iy)
        if (ix == 0 .and. iz == 0) e00 = e00 + w*dreal(cu*conjg(cu) + cv*conjg(cv) + cw*conjg(cw))
        call add(sum_new, w)
        if (.not. (ix == 0 .and. iz == 0)) call add(sum_old, w)
      end do
    end do
  end do
  do c = 1, 3
    call line_solve(KIND_DY, 0.0d0, V(:, :, :, c), rhs(:, :, :, 1))
    !$omp target update from(rhs)
    grad = 0.0d0
    do ix = nx0, nxN
      do iz = -nz, nz
        do iy = ny0, nyN
          w = 2.0d0*dyl(iy); if (ix == 0) w = dyl(iy)
          dq = rhs(iy, iz, ix, 1)
          grad = grad + w*dreal(dq*conjg(dq))
          if (.not. (ix == 0 .and. iz == 0)) sum_old(7) = sum_old(7) + w*dreal(dq*conjg(dq))
        end do
      end do
    end do
    sum_new(7) = sum_new(7) + grad
  end do
  sum_new(8) = e00
  call MPI_Allreduce(MPI_IN_PLACE, sum_new, 8, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_Allreduce(MPI_IN_PLACE, sum_old, 8, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  ! columns of the two files: energy diss uv vw ; uu vv ww uw ; E(0,0)
  call table(sum_new, ref); call table(sum_old, refold)
  ref(9) = 0.5d0*sum_new(8)

  passed = .true.
  if (has_terminal) then
    open (newunit=unit, file='statistics/Runtimedata', action='read'); read (unit, *) got_rt; close (unit)
    open (newunit=unit, file='statistics/variances_runtime.dat', action='read'); read (unit, *) got_var; close (unit)
    got_e00 = -1.0d0
    open (newunit=unit, file='exz_test.dat', action='read')
    read (unit, '(A)') line
    do
      read (unit, *, iostat=ios) ixx, izz, kx, kz, w
      if (ios /= 0) exit
      if (ixx == 0 .and. izz == 0) got_e00 = w
    end do
    close (unit)
    ! got:  energy diss uv vw  uu vv ww uw  E00
    errmax = 0.0d0
    call cmp('energy', got_rt(10), ref(1), refold(1))
    call cmp('diss  ', got_rt(11), ref(2), refold(2))
    call cmp('uv    ', got_rt(12), ref(3), refold(3))
    call cmp('vw    ', got_rt(13), ref(4), refold(4))
    call cmp('uu    ', got_var(2), ref(5), refold(5))
    call cmp('vv    ', got_var(3), ref(6), refold(6))
    call cmp('ww    ', got_var(4), ref(7), refold(7))
    call cmp('uw    ', got_var(5), ref(8), refold(8))
    call cmp('E(0,0)', got_e00, ref(9), 0.0d0)
    write (*, '(A,ES10.3,A,2ES10.3)') '  max relative error ', errmax, '   box means (u, w) removed: ', ubar, wbar
    passed = (errmax < 1.0d-10)
    d_old = abs(got_var(2) - refold(5))/abs(ref(5)) + abs(got_var(4) - refold(7))/abs(ref(7))
    if (d_old < 1.0d-3) then
      print *, 'the old operator gives the same uu, ww: the test does not see the change'
      passed = .false.
    end if
  end if
  call MPI_Bcast(passed, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
  call test_finish(passed)

contains

  ! this rank's mode / row summand w*(...) into  uu vv ww uv vw uw  eps_spectral  (index 7 also gets the k2 sum)
  subroutine add(a, wt)
    real(C_DOUBLE), intent(inout) :: a(8)
    real(C_DOUBLE), intent(in) :: wt
    a(1) = a(1) + wt*dreal(cu*conjg(cu)); a(2) = a(2) + wt*dreal(cv*conjg(cv)); a(3) = a(3) + wt*dreal(cw*conjg(cw))
    a(4) = a(4) + wt*dreal(cu*conjg(cv)); a(5) = a(5) + wt*dreal(cw*conjg(cv)); a(6) = a(6) + wt*dreal(cu*conjg(cw))
    a(7) = a(7) + wt*k2(iz, ix)*dreal(cu*conjg(cu) + cv*conjg(cv) + cw*conjg(cw))
  end subroutine add

  ! sums -> energy diss uv vw uu vv ww uw (box-mean units)
  subroutine table(a, r)
    real(C_DOUBLE), intent(in) :: a(8)
    real(C_DOUBLE), intent(out) :: r(9)
    r(1) = 0.5d0*(a(1) + a(2) + a(3))/ly
    r(2) = ni*a(7)/ly
    r(3) = a(4)/ly; r(4) = a(5)/ly
    r(5) = a(1)/ly; r(6) = a(2)/ly; r(7) = a(3)/ly; r(8) = a(6)/ly
    r(9) = 0.0d0
  end subroutine table

  subroutine cmp(name, got, expect, old)
    character(len=*), intent(in) :: name
    real(C_DOUBLE), intent(in) :: got, expect, old
    real(C_DOUBLE) :: e
    e = abs(got - expect)/max(abs(expect), 1.0d-300)
    errmax = max(errmax, e)
    write (*, '(2A,3(A,ES14.6))') '  ', name, '  live ', got, '   reference ', expect, '   old operator ', old
  end subroutine cmp

end program test_stats
