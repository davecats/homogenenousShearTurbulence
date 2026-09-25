! Everything the run knows: mesh, physical parameters, clock, rank layout,
! derived grid arrays and the solution fields.  Declarations only.
!
! This module must not `use` anything but the intrinsic iso_c_binding.
! nvfortran fails to resolve `declare target` variables inside downstream
! target regions when the module declaring them also `use`s other modules
! (see channel/src/physics/README.md).  Anything that needs a `use` goes in
! another module and fills these variables in.
!
! Axes: x streamwise (Fourier, alfa0), z spanwise (Fourier, beta0),
!       y the compact-finite-difference direction carrying the shear U = S*y.
! Field layout: V(iy, iz, ix, i) with i = 1:u, 2:v, 3:w ; iy = -2:ny+1
!       (unique points 0..ny-1, four shear-periodic ghost rows), iz = -nz:nz,
!       ix = nx0:nxN (this rank's x slab).
module hst_params

  use, intrinsic :: iso_c_binding

  implicit none

  !---------------------------------------------------------------- mesh ----
  integer(C_INT), save :: nx, ny, nz          ! modes in x, points in y, modes in z
  integer(C_INT), save :: nxd, nzd            ! 3/2-padded transform sizes
  real(C_DOUBLE), save :: alfa0, beta0        ! fundamental wavenumbers
  real(C_DOUBLE), save :: lx, ly, lz          ! box; ly is an input, lx, lz derived
  real(C_DOUBLE), save :: ystretch            ! tanh clustering of y at mid-box (0 = uniform)
  integer(C_INT), save :: line_chunk          ! x columns per line-solver batch (0 = all, hst_linsolve)
  !$omp declare target(ny)

  !------------------------------------------------------------- physics ----
  real(C_DOUBLE), save :: re, ni              ! Reynolds number and 1/re
  real(C_DOUBLE), save :: S                   ! mean shear dU/dy
  logical, save :: linear                     ! .true.: nonlinear terms off
  logical, save :: exact_shift                ! .true.: advect the unweighted quantities (FINDINGS.md)
  ! unsteady spanwise mean shear dW/dy = S2(t) = s2_amplitude * sin(2 pi (t - s2_start)/s2_period)
  ! for t >= s2_start (constant s2_amplitude when s2_period = 0); zero amplitude switches it off
  real(C_DOUBLE), save :: s2_amplitude, s2_period, s2_start
  ! Stokes layer: oscillating spanwise mean profile at mid-box (hst_stokes)
  real(C_DOUBLE), save :: sl_amplitude, sl_period, sl_delta, sl_start
  logical, save :: sl_bodyforce, sl_ramp
  !$omp declare target(ni, S)

  !------------------------------------------------------- clock and I/O ----
  real(C_DOUBLE), save :: deltat, dt_fixed, cflmax, cfl = 0.0d0   ! dt_fixed: deck value, 0 = from cflmax
  real(C_DOUBLE), save :: t_max, dt_stat, dt_field, dt_save
  real(C_DOUBLE), save :: time, time0 = 0.0d0
  integer(C_SIZE_T), save :: nstep, istep = 0, ifield = 0
  logical, save :: time_from_restart

  !------------------------------------------------------ initial field ----
  real(C_DOUBLE), save :: amplitude, kpeak
  integer(C_INT), save :: seed

  !-------------------------------------------------------- rank layout ----
  ! nproc ranks form a 1-D x-z pencil decomposition (npxz = nproc, npy = 1).
  ! The y range ny0:nyN is kept as in the channel code so that a y
  ! decomposition can be added later; today ny0 = 0 and nyN = ny-1.
  integer(C_INT), save :: nproc, iproc, npxz, ipxz
  integer(C_INT), save :: nx0, nxN, nxB, nz0, nzN, nzB, ny0, nyN
  logical, save :: has_terminal, has_average
  !$omp declare target(nx0, nxN, nxB, nz0, nzN, nzB, ny0, nyN)

  !--------------------------------------------------- derived grid data ----
  real(C_DOUBLE), allocatable, save :: y(:)             ! y(-2:ny+1)
  real(C_DOUBLE), allocatable, save :: dyl(:)           ! dyl(0:ny-1): row spacing 0.5 (y(iy+1) - y(iy-1))
  real(C_DOUBLE), allocatable, save :: fy(:)            ! fy(-2:ny+1): body force on the mean w (Stokes layer)
  integer(C_INT), allocatable, save :: inlayer(:)       ! inlayer(0:ny-1): 1 inside the Stokes layer
  real(C_DOUBLE), allocatable, save :: der(:, :, :)     ! der(0:ny-1, 0:3, -2:2), see hst_derivatives
  real(C_DOUBLE), allocatable, save :: k2(:, :)         ! alfa^2+beta^2 (iz, ix)
  complex(C_DOUBLE_COMPLEX), allocatable, save :: ialfa(:), ibeta(:)
  integer(C_INT), allocatable, save :: izd(:)           ! z mode -> padded index
  real(C_DOUBLE), save :: dx, dy, dz, factor            ! CFL spacings, 1/(2 nxd nzd)

  !---------------------------------------------------- solution fields ----
  ! V: the velocity, with ghost rows.  oldrhs: the explicit terms of the
  ! previous substep, i = 1:eta, 2:d2v.  rhs: the right-hand sides of the two
  ! equations inside a substep (buildrhs_prepare .. linsolve, same i); free
  ! between substeps, where it is scratch for outstats (i = 1) and for the
  ! pressure (i = 1 right-hand side, i = 2 result).
  complex(C_DOUBLE_COMPLEX), allocatable, target, save :: V(:, :, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable, save :: oldrhs(:, :, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable, save :: rhs(:, :, :, :)

  !--------------------------------------- Runge-Kutta (Rai-Moin) table ----
  ! One column per substep: (1) unknown/deltat, (2) new explicit, (3) old
  ! explicit.  Substep lengths are 2/RK_rai(1, i)*deltat = 8/15, 2/15, 1/3.
  real(C_DOUBLE), dimension(3, 3), save :: RK_rai = reshape( &
      (/120.0d0/32.0d0, 2.0d0, 0.0d0, &
        120.0d0/8.0d0, 50.0d0/8.0d0, 34.0d0/8.0d0, &
        120.0d0/20.0d0, 90.0d0/20.0d0, 50.0d0/20.0d0/), shape=(/3, 3/))

end module hst_params
