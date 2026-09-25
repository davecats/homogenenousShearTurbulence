! Kelvin-mode test of the linearised solver (WP2, extended to S2).
!
!   mpirun -np 1 build-cpu/test_kelvin [hst.in]     (linear = .true. is forced)
!
! A single Fourier mode in the uniform shear U = S y (+ W = S2(t) y) is an
! exact solution of the linearised equations (Kelvin 1887; Moffatt 1967):
! its wavevector tilts,
!   ky(t) = ky0 - S kx t - kz gamma_y(t),      gamma_y = int S2 dt,
! the vertical velocity follows  v(t) = v0 k0^2/k(t)^2  and the vertical
! vorticity  eta(t) = eta0 + int_0^t (S2 i kx - S i kz) v dt' , both times
! the viscous factor  exp(-ni int_0^t k^2 dt').  With the time-dependent
! shear-periodic condition the mode stays a single mode in x, z with y
! dependence exp(i ky(t) y), so this exercises the wrap phases, the exact
! mean-shear advection, the tilting terms and the implicit solves at once.
!
! v is checked for any S2(t) (the k^2 integral is done numerically); eta
! only for constant S2 (s2_period = 0), where its integral is closed-form.
! Discretisation errors: O(deltat^3) in time, and with the default
! advection the (1/6) dy^2 (ky0^2 - ky(t)^2) term of PLAN.md 8; with
! exact_shift only the sixth-order stencil error remains.  Passes at 1e-3
! relative (tests/decks/kelvin*.in).
program test_kelvin

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

  implicit none

  integer :: ierr, iy, m, n, iq
  real(C_DOUBLE), parameter :: PI = 3.141592653589793d0
  complex(C_DOUBLE_COMPLEX), parameter :: I = (0.0d0, 1.0d0)
  integer(C_INT), parameter :: ix = 1, iz = 1, mky = 1
  complex(C_DOUBLE_COMPLEX), parameter :: v0 = (1.0d0, 0.0d0), eta0 = (0.3d0, -0.2d0)
  real(C_DOUBLE) :: kx, kz, ky0, kh, kh2, k02, kyt, k2t, visc, err_v, err_eta, ratio, t_end
  real(C_DOUBLE) :: c, tq, k2int, s2c
  complex(C_DOUBLE_COMPLEX) :: vex, etaex, u0, w0, ph
  logical :: check_eta

  call test_start()
  linear = .true.
  time = 0.0d0
  if (nproc /= 1) error stop 'test_kelvin runs on one rank'
  call test_setup()

  ! The mode at t = 0: v = v0 e^{i ky0 y}, u and w from continuity and eta0.
  kx = alfa0*ix; kz = beta0*iz; ky0 = 2.0d0*PI*mky/ly
  kh2 = kx*kx + kz*kz; kh = sqrt(kh2)
  k02 = kh2 + ky0*ky0
  u0 = (I*kx*(I*ky0*v0) - I*kz*eta0)/kh2
  w0 = (I*kz*(I*ky0*v0) + I*kx*eta0)/kh2
  V = 0
  do iy = -2, ny + 1
    ph = exp(I*ky0*y(iy))
    V(iy, iz, ix, 1) = u0*ph
    V(iy, iz, ix, 2) = v0*ph
    V(iy, iz, ix, 3) = w0*ph
  end do
  !$omp target update to(V)

  ! Advance with a fixed step to t_end.
  if (dt_fixed <= 0.0d0) dt_fixed = 1.0d-3
  deltat = dt_fixed
  t_end = t_max
  n = nint(t_end/deltat)
  do m = 1, n
    call timestep()
  end do
  !$omp target update from(V)

  ! Closed form at time: ky(t), the viscous factor (k^2 integrated by the
  ! trapezoidal rule on a fine grid), v, and eta for constant S2.
  kyt = ky0 - S*kx*time - kz*gamma_y_of(time)
  k2t = kh2 + kyt*kyt
  k2int = 0.0d0
  do iq = 0, 20000
    tq = time*iq/20000.0d0
    c = kh2 + (ky0 - S*kx*tq - kz*gamma_y_of(tq))**2
    if (iq == 0 .or. iq == 20000) c = 0.5d0*c
    k2int = k2int + c
  end do
  k2int = k2int*time/20000.0d0
  visc = exp(-ni*k2int)
  vex = v0*k02/k2t*visc
  check_eta = (s2_period == 0.0d0 .and. s2_start <= 0.0d0)
  s2c = s2_of(0.0d0)
  c = S*kx + s2c*kz                        ! d ky/dt = -c for constant S2
  if (abs(c) > 0.0d0) then
    etaex = (eta0 + (s2c*I*kx - S*I*kz)*v0*k02/(kh*c)*(atan(ky0/kh) - atan(kyt/kh)))*visc
  else
    etaex = (eta0 + (s2c*I*kx - S*I*kz)*v0*time)*visc
  end if
  err_v = 0; err_eta = 0
  do iy = 0, ny - 1
    ph = exp(I*kyt*y(iy))
    err_v = max(err_v, abs(V(iy, iz, ix, 2) - vex*ph))
    err_eta = max(err_eta, abs(I*kz*V(iy, iz, ix, 1) - I*kx*V(iy, iz, ix, 3) - etaex*ph))
  end do
  ratio = maxval(abs(V(:, iz, ix, 2)))/abs(v0)
  if (has_terminal) then
    write (*, '(A,F8.4,A,I0,A,ES9.2,A,F6.3,A,F6.3)') '   Kelvin mode to t = ', time, ' in ', n, ' steps of ', deltat, &
      '   S = ', S, '  S2(0) = ', s2c
    write (*, '(A,F8.4,A,F8.4,A,F8.4)') '   ky: ', ky0, ' -> ', kyt, '   |v|/|v0| = ', ratio
    write (*, '(A,ES10.2,A,ES10.2)') '   max error in v:   ', err_v, '   relative ', err_v/abs(vex)
    if (check_eta) then
      write (*, '(A,ES10.2,A,ES10.2)') '   max error in eta: ', err_eta, '   relative ', err_eta/abs(etaex)
    else
      write (*, '(A)') '   (eta not checked: no closed form for time-dependent S2)'
    end if
  end if
  call test_finish(err_v/abs(vex) <= 1.0d-3 .and. (.not. check_eta .or. err_eta/abs(etaex) <= 1.0d-3))

end program test_kelvin
