! Stokes-layer test (section 10 wish).
!
!   mpirun -np N build-cpu/test_stokes tests/decks/stokes.in
!
! With no fluctuations the mean spanwise profile obeys  dW/dt = nu d2W/dy2 + f
! and f was built so that the analytic Stokes layer W(y, t) is its solution.
! The profile is prescribed during the first two periods after sl_start;
! the deck runs one period longer, so the last period is driven by the body
! force and the viscous term alone.  The (0,0) mode of w at the end must
! equal W(y, t) to the discretisation error (CN in time, sixth-order
! stencils on the stretched grid); everything else must stay zero.
program test_stokes

  use, intrinsic :: iso_c_binding
  use mpi_f08
  use hst_params
  use test_common
  use hst_input
  use hst_mpi
  use hst_fft
  use hst_setup
  use hst_derivatives
  use hst_linsolve
  use hst_equations
  use hst_stokes
  use hst_io

  implicit none

  integer :: ierr, iy, m, n
  real(C_DOUBLE) :: err, wmax, other, err_g, other_g, wmax_g, wex

  call test_start()
  if (sl_amplitude == 0.0d0) error stop 'test_stokes needs a deck with sl_amplitude /= 0'
  time = 0.0d0
  call test_setup()
  V = 0
  !$omp target update to(V)
  call stokes_setup()
  call stokes_apply()
  do m = 1, 3
    call fill_ghosts(m)
  end do

  if (dt_fixed <= 0.0d0) dt_fixed = 5.0d-3
  deltat = dt_fixed
  n = nint(t_max/deltat)
  do m = 1, n
    call timestep()
  end do
  !$omp target update from(V)

  err = 0; wmax = 0; other = 0
  if (has_average) then
    do iy = 0, ny - 1
      wex = stokes_profile(y(iy), time)
      err = max(err, abs(V(iy, 0, 0, 3) - wex))
      wmax = max(wmax, abs(wex))
      other = max(other, abs(dimag(V(iy, 0, 0, 3))), abs(V(iy, 0, 0, 1)), abs(V(iy, 0, 0, 2)))
    end do
    V(:, 0, 0, 3) = 0
  end if
  other = max(other, maxval(abs(V(0:ny - 1, :, :, :))))
  call MPI_Allreduce(err, err_g, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  call MPI_Allreduce(other, other_g, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  call MPI_Allreduce(wmax, wmax_g, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  wmax = wmax_g
  if (has_terminal) then
    write (*, '(A,F7.3,A,I0,A,F6.3,A,F6.3)') '   Stokes layer to t = ', time, ' in ', n, ' steps; delta = ', sl_delta, &
      '  min dy = ', minval(dyl)
    write (*, '(A,ES10.2,A,ES10.2,A,ES10.2)') '   max error of the mean profile ', err_g, '  (max |W| = ', wmax, &
      ')  relative ', err_g/wmax
    write (*, '(A,ES10.2)') '   everything else (must be zero): ', other_g
  end if
  call test_finish(err_g/wmax <= 1.0d-4 .and. other_g <= 1.0d-12)

end program test_stokes
