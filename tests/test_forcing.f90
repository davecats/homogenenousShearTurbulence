! Nonlinear forcing of the vertical-vorticity equation (WP5).
!
!   mpirun -np 1 build-cpu/test_forcing [hst.in]     (S = 0, inviscid, forced)
!
! For  u = sin(ky),  v = 0,  w = sin(kx)  (divergence-free), the nonlinear
! term is  N_u = 0,  N_v = 0,  N_w = -k sin(ky) cos(kx), so
!   d eta/dt = dN_u/dz - dN_w/dx = -k^2 sin(ky) sin(kx),   d(lap v)/dt = 0.
! In modes: sin(kx) has coefficient -i/2 at ix = mx, so after one inviscid
! step of length dt the eta mode (mx, 0) must change by
!   dt * (-k^2) * (-i/2) * sin(ky) = i k^2 dt/2 sin(ky)
! to O(dt^2), and v must stay zero.  u is a mean profile (the (0,0) mode),
! so this also exercises the packed mean-mode path.  k = 2 pi: mx = 3
! with alfa0 = 2 pi/3, my = 2 with ly = 2.
program test_forcing

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
  use hst_transforms
  use hst_equations

  implicit none

  integer :: ierr, iy, m
  real(C_DOUBLE), parameter :: PI = 3.141592653589793d0
  complex(C_DOUBLE_COMPLEX), parameter :: I = (0.0d0, 1.0d0)
  integer(C_INT), parameter :: mx = 3, my = 2
  real(C_DOUBLE) :: k, ky, err, errv, ref
  complex(C_DOUBLE_COMPLEX) :: eta0, eta1, expected

  call test_start()
  S = 0.0d0; linear = .false.; ni = 1.0d-12; time = 0.0d0
  if (nproc /= 1) error stop 'test_forcing runs on one rank'
  if (abs(mx*alfa0 - 2.0d0*PI*my/ly) > 1.0d-12) error stop 'test_forcing needs alfa0 = 2pi/3, ly = 2'
  call test_setup()

  k = mx*alfa0; ky = k
  V = 0
  do iy = -2, ny + 1
    V(iy, 0, 0, 1) = sin(ky*y(iy))            ! u mean profile
    V(iy, 0, mx, 3) = -0.5d0*I                ! w = sin(kx)
  end do
  !$omp target update to(V)
  do m = 1, 3
    call fill_ghosts(m)
  end do
  deltat = 1.0d-4
  eta0 = ibeta(0)*V(0, 0, mx, 1) - ialfa(mx)*V(0, 0, mx, 3)
  call timestep()
  !$omp target update from(V)

  err = 0.0d0; errv = 0.0d0; ref = 0.0d0
  do iy = 0, ny - 1
    eta1 = ibeta(0)*V(iy, 0, mx, 1) - ialfa(mx)*V(iy, 0, mx, 3)
    expected = eta0 + I*k*k*deltat/2.0d0*sin(ky*y(iy))
    err = max(err, abs(eta1 - expected))
    ref = max(ref, abs(expected - eta0))
    errv = max(errv, abs(V(iy, 0, mx, 2)))
  end do
  if (has_terminal) then
    write (*, '(A,ES10.2,A,ES10.2,A,ES10.2)') '   eta forcing: change ', ref, '  error ', err, '  relative ', err/ref
    write (*, '(A,ES10.2)') '   v (must stay zero): ', errv
  end if
  call test_finish(err/ref <= 1.0d-2 .and. errv <= 1.0d-8)

end program test_forcing
