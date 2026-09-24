! Kelvin-mode test of the linearised solver (WP2).
!
!   mpirun -np 1 build-cpu/test_kelvin [hst.in]     (linear = .true. is forced)
!
! A single Fourier mode in uniform shear U = S y is an exact solution of the
! linearised equations (Kelvin 1887; Moffatt 1967): its wavevector tilts,
!   ky(t) = ky0 - S kx t,
! the vertical velocity follows  v(t) = v0 k0^2/k(t)^2  and the vertical
! vorticity  eta(t) = eta0 - i S kz int_0^t v dt' , both times the viscous
! factor  exp(-ni int_0^t k^2 dt').  With the time-dependent shear-periodic
! condition the mode stays a single mode in x, z with y dependence
! exp(i ky(t) y), so this exercises the wrap phase, the exact mean-shear
! advection, the tilting term and the implicit solves at once.
!
! Discretisation errors: O(deltat^3) in time, and the compact scheme's
! modified wavenumber in y, O((ky dy)^4).  The test passes when the field
! agrees with the closed form to 1e-3 relative, and prints the error so that
! halving dy can be seen to divide it by ~16.
program test_kelvin

  use, intrinsic :: iso_c_binding
  use mpi_f08
  use hst_params
  use hst_input
  use hst_mpi
  use hst_fft
  use hst_setup
  use hst_derivatives
  use hst_linsolve
  use hst_equations
#ifdef HAVE_CUDA
  use omp_lib
#endif

  implicit none

  character(len=256) :: deck
  integer :: ierr, iy, m, n
  real(C_DOUBLE), parameter :: PI = 3.141592653589793d0
  complex(C_DOUBLE_COMPLEX), parameter :: I = (0.0d0, 1.0d0)
  integer(C_INT), parameter :: ix = 1, iz = 1, mky = 1
  complex(C_DOUBLE_COMPLEX), parameter :: v0 = (1.0d0, 0.0d0), eta0 = (0.3d0, -0.2d0)
  real(C_DOUBLE) :: kx, kz, ky0, kh2, k02, kyt, k2t, visc, err_v, err_eta, ratio, t_end
  complex(C_DOUBLE_COMPLEX) :: vex, etaex, u0, w0, ph

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, iproc, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
  has_terminal = (iproc == 0)
#ifdef HAVE_CUDA
  call omp_set_default_device(mod(iproc, omp_get_num_devices()))
#endif
  deck = 'hst.in'
  if (command_argument_count() >= 1) call get_command_argument(1, deck)
  call read_input(trim(deck))
  linear = .true.
  time = 0.0d0
  if (nproc /= 1) error stop 'test_kelvin runs on one rank'
  call setup_decomposition()
  call allocate_fields()
  call init_fft()
  call init_linsolve()
  call setup_derivatives()

  ! The mode at t = 0: v = v0 e^{i ky0 y}, u and w from continuity and eta0.
  kx = alfa0*ix; kz = beta0*iz; ky0 = 2.0d0*PI*mky/ly
  kh2 = kx*kx + kz*kz
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

  ! Closed form at time
  kyt = ky0 - S*kx*time
  k2t = kh2 + kyt*kyt
  visc = exp(-ni*(kh2*time + ky0*ky0*time - ky0*S*kx*time*time + S*S*kx*kx*time**3/3.0d0))
  vex = v0*k02/k2t*visc
  if (abs(S*kx) > 0.0d0) then
    etaex = (eta0 - I*S*kz*v0*k02/(sqrt(kh2)*S*kx)*(atan(ky0/sqrt(kh2)) - atan(kyt/sqrt(kh2))))*visc
  else
    etaex = (eta0 - I*S*kz*v0*time)*visc
  end if
  err_v = 0; err_eta = 0
  do iy = 0, ny - 1
    ph = exp(I*kyt*y(iy))
    err_v = max(err_v, abs(V(iy, iz, ix, 2) - vex*ph))
    err_eta = max(err_eta, abs(I*kz*V(iy, iz, ix, 1) - I*kx*V(iy, iz, ix, 3) - etaex*ph))
  end do
  ratio = maxval(abs(V(:, iz, ix, 2)))/abs(v0)
  if (has_terminal) then
    write (*, '(A,F8.4,A,I0,A,ES9.2)') '   Kelvin mode to t = ', time, ' in ', n, ' steps of ', deltat
    write (*, '(A,F8.4,A,F8.4,A,F8.4)') '   ky: ', ky0, ' -> ', kyt, '   |v|/|v0| = ', ratio
    write (*, '(A,ES10.2,A,ES10.2)') '   max error in v:   ', err_v, '   relative ', err_v/abs(vex)
    write (*, '(A,ES10.2,A,ES10.2)') '   max error in eta: ', err_eta, '   relative ', err_eta/abs(etaex)
  end if
  call free_linsolve(); call free_fft(); call free_fields(); call free_mpi()
  call MPI_Finalize(ierr)
  if (err_v/abs(vex) > 1.0d-3 .or. err_eta/abs(etaex) > 1.0d-3) then
    print *, 'FAILED'
    error stop 1
  end if
  if (has_terminal) print *, 'PASSED'

end program test_kelvin
