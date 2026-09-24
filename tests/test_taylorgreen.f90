! Taylor-Green test of the full nonlinear solver (WP5).
!
!   mpirun -np 1 build-cpu/test_taylorgreen [hst.in]     (S = 0 is forced)
!
! A two-dimensional Taylor-Green vortex with equal wavenumbers is an exact
! solution of the Navier-Stokes equations: its nonlinear term is a pure
! gradient (absorbed by the pressure), so every mode decays as
! exp(-ni (kx^2 + ky^2) t) with the nonlinear terms switched ON.  Two
! orientations are run, each for n fixed steps from time zero:
!   x-y plane:  u = cos(kx) sin(ky),  v = -sin(kx) cos(ky)      (uu, uv, vv, alfa terms)
!   z-y plane:  w = cos(kz) sin(ky),  v = -sin(kz) cos(ky)      (ww, vw, vv, beta terms)
! with k = 2 pi (ix = 3 with alfa0 = 2 pi/3; iz = 1 with beta0 = 2 pi; ky
! with my = 2 for ly = 2).  Passes when the field agrees with the exact
! decay to 1e-4 relative (6.7e-6 measured, time-integration error); a wrong
! nonlinear term shows up at O(1).
program test_taylorgreen

  use, intrinsic :: iso_c_binding
  use mpi_f08
  use hst_params
  use hst_input
  use hst_mpi
  use hst_fft
  use hst_setup
  use hst_derivatives
  use hst_linsolve
  use hst_transforms
  use hst_equations
#ifdef HAVE_CUDA
  use omp_lib
#endif

  implicit none

  character(len=256) :: deck
  integer :: ierr, iy, ix, iz, m, n, orient
  real(C_DOUBLE), parameter :: PI = 3.141592653589793d0
  complex(C_DOUBLE_COMPLEX), parameter :: I = (0.0d0, 1.0d0)
  integer(C_INT), parameter :: mx = 3, mz = 1, my = 2
  real(C_DOUBLE) :: k, ky, decay, err, err_global, worst
  complex(C_DOUBLE_COMPLEX) :: ex

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
  S = 0.0d0; linear = .false.
  if (nx < 2*mx) error stop 'test_taylorgreen needs nx >= 6'
  if (abs(mx*alfa0 - 2.0d0*PI*my/ly) > 1.0d-12 .or. abs(mz*beta0 - 2.0d0*PI*my/ly) > 1.0d-12) &
    error stop 'test_taylorgreen needs alfa0 = 2pi/3, beta0 = 2pi, ly = 2'
  call setup_decomposition()
  call allocate_fields()
  call init_fft()
  call init_linsolve()
  call setup_derivatives()
  if (dt_fixed <= 0.0d0) dt_fixed = 1.0d-3
  n = 20
  worst = 0.0d0

  do orient = 1, 2
    k = 2.0d0*PI; ky = k
    time = 0.0d0; deltat = dt_fixed; oldrhs = 0
    !$omp target update to(oldrhs)
    V = 0
    if (orient == 1) then          ! cos(kx) sin(ky), -sin(kx) cos(ky): mode ix = mx, iz = 0
      if (nx0 <= mx .and. mx <= nxN) then
        do iy = -2, ny + 1
          V(iy, 0, mx, 1) = 0.5d0*sin(ky*y(iy))
          V(iy, 0, mx, 2) = 0.5d0*I*cos(ky*y(iy))
        end do
      end if
    else                           ! cos(kz) sin(ky), -sin(kz) cos(ky): modes ix = 0, iz = +-mz
      if (nx0 == 0) then
        do iy = -2, ny + 1
          V(iy, mz, 0, 3) = 0.5d0*sin(ky*y(iy));   V(iy, -mz, 0, 3) = 0.5d0*sin(ky*y(iy))
          V(iy, mz, 0, 2) = 0.5d0*I*cos(ky*y(iy)); V(iy, -mz, 0, 2) = -0.5d0*I*cos(ky*y(iy))
        end do
      end if
    end if
    !$omp target update to(V)
    do m = 1, 3
      call fill_ghosts(m)
    end do
    do m = 1, n
      call timestep()
    end do
    !$omp target update from(V)

    decay = exp(-ni*2.0d0*k*k*time)
    err = 0.0d0
    do ix = nx0, nxN
      do iz = -nz, nz
        do iy = 0, ny - 1
          do m = 1, 3
            ex = 0.0d0
            if (orient == 1 .and. ix == mx .and. iz == 0) then
              if (m == 1) ex = 0.5d0*sin(ky*y(iy))*decay
              if (m == 2) ex = 0.5d0*I*cos(ky*y(iy))*decay
            else if (orient == 2 .and. ix == 0 .and. abs(iz) == mz) then
              if (m == 3) ex = 0.5d0*sin(ky*y(iy))*decay
              if (m == 2) ex = sign(1, iz)*0.5d0*I*cos(ky*y(iy))*decay
            end if
            err = max(err, abs(V(iy, iz, ix, m) - ex))
          end do
        end do
      end do
    end do
    call MPI_Allreduce(err, err_global, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    worst = max(worst, err_global/(0.5d0*decay))
    if (has_terminal) write (*, '(A,I1,A,I0,A,ES10.2,A,ES10.2)') '   Taylor-Green orientation ', orient, &
      ': after ', n, ' steps max error ', err_global, '  relative ', err_global/(0.5d0*decay)
  end do

  call free_linsolve(); call free_fft(); call free_fields(); call free_mpi()
  call MPI_Finalize(ierr)
  if (worst > 1.0d-4) then
    print *, 'FAILED'
    error stop 1
  end if
  if (has_terminal) print *, 'PASSED'

end program test_taylorgreen
