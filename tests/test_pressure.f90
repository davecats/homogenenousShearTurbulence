! Pressure test (WP3b): Taylor-Green vortex in the x-y plane.
!
!   mpirun -np 1 build-cpu/test_pressure [hst.in]     (S = 0 is forced)
!
! For  u = cos(kx) sin(ky),  v = -sin(kx) cos(ky),  w = 0  with equal
! wavenumbers k, the pressure is  p = -(cos 2kx + cos 2ky)/4:  in modes,
! p(ix = 2*mx, iz = 0) = -1/8 for all y and p(0, 0)(y) = -cos(2ky)/4 (the
! latter is the singular mean mode, fixed by the zero-mean gauge).  The
! test passes at 1e-4 (1.9e-5 at ny = 64, fourth-order convergence).  It
! uses k = 2 pi: mx = 3 (kx = 3 alfa0 with the default alfa0 = 2 pi/3)
! and my = 2 (ky = 2 pi my/ly with ly = 2), so the deck must have nx >= 6.
! The x-y plane exercises the D1 uv and D2 vv terms and the mean mode.
program test_pressure

  use, intrinsic :: iso_c_binding
  use mpi_f08
  use hst_params
  use hst_input
  use hst_mpi
  use hst_fft
  use hst_setup
  use hst_derivatives
  use hst_linsolve
  use hst_pressure
#ifdef HAVE_CUDA
  use omp_lib
#endif

  implicit none

  character(len=256) :: deck
  integer :: ierr, iy, ix, iz, m
  real(C_DOUBLE), parameter :: PI = 3.141592653589793d0
  complex(C_DOUBLE_COMPLEX), parameter :: I = (0.0d0, 1.0d0)
  integer(C_INT), parameter :: mx = 3, my = 2
  real(C_DOUBLE) :: k, ky, err, err_other, pex, err_global
  complex(C_DOUBLE_COMPLEX) :: ep, em

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
  S = 0.0d0; time = 0.0d0
  if (nx < 2*mx) error stop 'test_pressure needs nx >= 6'
  if (abs(mx*alfa0 - 2.0d0*PI*my/ly) > 1.0d-12) error stop 'test_pressure needs kx = ky: alfa0 = 2pi/3, ly = 2'
  call setup_decomposition()
  call allocate_fields()
  call init_fft()
  call init_linsolve()
  call setup_derivatives()

  ! cos(kx) sin(ky) = (e^{ikx} + e^{-ikx})/2 * sin(ky): stored mode ix = mx
  ! with coefficient 1/2 (the -mx partner is implied by Hermitian symmetry);
  ! -sin(kx) cos(ky) = -(e^{ikx} - e^{-ikx})/(2i) cos(ky): coefficient -1/(2i) = i/2.
  k = mx*alfa0; ky = k
  V = 0
  if (nx0 <= mx .and. mx <= nxN) then
    do iy = -2, ny + 1
      V(iy, 0, mx, 1) = 0.5d0*sin(ky*y(iy))
      V(iy, 0, mx, 2) = 0.5d0*I*cos(ky*y(iy))
    end do
  end if
  !$omp target update to(V)
  do m = 1, 3
    call fill_ghosts(m)
  end do

  call compute_pressure(rhs(:, :, :, 2))
  !$omp target update from(rhs)

  err = 0; err_other = 0
  do ix = nx0, nxN
    do iz = -nz, nz
      do iy = 0, ny - 1
        if (ix == 2*mx .and. iz == 0) then
          err = max(err, abs(rhs(iy, iz, ix, 2) - dcmplx(-0.125d0, 0.0d0)))
        else if (ix == 0 .and. iz == 0) then
          pex = -0.25d0*cos(2.0d0*ky*y(iy))
          err = max(err, abs(rhs(iy, iz, ix, 2) - pex))
        else
          err_other = max(err_other, abs(rhs(iy, iz, ix, 2)))
        end if
      end do
    end do
  end do
  call MPI_Allreduce(err, err_global, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  err = err_global
  call MPI_Allreduce(err_other, err_global, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  err_other = err_global
  if (has_terminal) then
    write (*, '(A,ES10.2)') '   Taylor-Green pressure: max error in the two nonzero modes ', err
    write (*, '(A,ES10.2)') '                          max magnitude of all other modes    ', err_other
  end if
  call free_linsolve(); call free_fft(); call free_fields(); call free_mpi()
  call MPI_Finalize(ierr)
  if (err > 1.0d-4 .or. err_other > 1.0d-10) then
    print *, 'FAILED'
    error stop 1
  end if
  if (has_terminal) print *, 'PASSED'

end program test_pressure
