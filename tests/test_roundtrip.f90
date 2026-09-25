! Round-trip test of the transform and I/O layers (WP1).
!
!   mpirun -np N build-cpu/test_roundtrip [hst.in]
!
! 1. A random field goes spectral -> physical -> spectral through the same
!    path the solver uses (assemble, IFT, transpose, RFT; then HFT,
!    transpose, FFT) and must come back unchanged, on every rank count.
! 2. The field is written with restart_write, read back with restart_read,
!    and must match.
! Prints the largest error of each; exits with status 1 if any exceeds 1e-12
! relative to the field's largest value.
program test_roundtrip

  use, intrinsic :: iso_c_binding
  use mpi_f08
  use hst_params
  use test_common
  use hst_input
  use hst_mpi
  use hst_fft
  use hst_setup
  use hst_transforms
  use hst_initial
  use hst_io

  implicit none

  integer :: ierr, m, i, j, k, y_first, y_last
  complex(C_DOUBLE_COMPLEX), allocatable :: V0(:, :, :, :), W(:, :, :)
  real(C_DOUBLE) :: err, err_global, vmax, tol, worst

  call test_start()
  call test_setup()

  ! A random field with plain periodic ghost rows (gamma = 0 at time zero).
  call generate_initial_field()
  V(-2:-1, :, :, :) = V(ny - 2:ny - 1, :, :, :)
  V(ny:ny + 1, :, :, :) = V(0:1, :, :, :)
  allocate (V0, source=V)
  allocate (W(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)); W = 0
  !$omp target enter data map(to: W)
  !$omp target update to(V)
  vmax = maxval(abs(V0))
  call MPI_Allreduce(MPI_IN_PLACE, vmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  tol = 1.0d-12*vmax
  worst = 0.0d0

  ! 1. spectral -> physical -> spectral
  call transform_to_physical()
  call compute_cfl()
  y_first = ny0 - 2
  y_last = nyN + 2
  do m = 1, 3
    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(rVVdx, products, nxd, nzB, y_first, y_last, factor, m) private(i, j, k)
    do i = y_first, y_last
      do j = 1, nzB
        do k = 1, 2*nxd
          products(k, j, i) = rVVdx(k, j, i, m)*factor
        end do
      end do
    end do
    call products_to_spectral()
    call vvdz_to_field(W)
    !$omp target update from(W)
    err = maxval(abs(W - V0(:, :, :, m)))
    call MPI_Allreduce(err, err_global, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    worst = max(worst, err_global)
    if (has_terminal) write (*, '(A,I1,A,ES10.2,A,ES10.2)') '   round trip component ', m, ': max error ', err_global, &
      '  (max |V| = ', vmax, ')'
  end do

  ! 2. write, read back
  call restart_write('test_roundtrip.field')
  V = 0
  call restart_read('test_roundtrip.field')
  err = maxval(abs(V(0:ny - 1, :, :, :) - V0(0:ny - 1, :, :, :)))
  call MPI_Allreduce(err, err_global, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  worst = max(worst, err_global)
  if (has_terminal) write (*, '(A,ES10.2)') '   restart write/read: max error ', err_global
  if (has_terminal) write (*, '(A,F10.4)') '   rank-0 cfl estimate ', cfl

  !$omp target exit data map(delete: W)
  call test_finish(worst <= tol)

end program test_roundtrip
