! Unit test of the cyclic pentadiagonal line solver (WP2).
!
!   mpirun -np 1 build-cpu/test_linsolve [hst.in]
!
! Random complex band matrices with the wrap-around structure of the
! shear-periodic stencil (corner entries carrying a phase), random exact
! solutions, right-hand sides formed by the band product; the solver must
! recover the solutions to round-off.  Runs on the device in a GPU build.
program test_linsolve

  use, intrinsic :: iso_c_binding
  use mpi_f08
  use hst_params
  use hst_input
  use hst_mpi
  use hst_linsolve
  use hst_initial, only: uniform_from_key
#ifdef HAVE_CUDA
  use omp_lib
#endif

  implicit none

  character(len=256) :: deck
  integer :: ierr, il, iy, j, c, nl
  complex(C_DOUBLE_COMPLEX), allocatable :: xtrue(:, :), Aref(:, :, :)
  complex(C_DOUBLE_COMPLEX) :: ph, b
  real(C_DOUBLE) :: err, xmax, r1, r2

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
  call setup_decomposition()
  call init_linsolve()

  nl = min(nlines_max, 37)
  allocate (xtrue(nl, 0:ny - 1), Aref(nl, 0:ny - 1, -2:2))
  do il = 1, nl
    ph = exp(dcmplx(0.0d0, 6.283185307179586d0*uniform_from_key(7, 9, il, 0, 0)))
    do iy = 0, ny - 1
      do j = -2, 2
        r1 = uniform_from_key(1, j + 3, iy, il, 0) - 0.5d0
        r2 = uniform_from_key(2, j + 3, iy, il, 0) - 0.5d0
        Aref(il, iy, j) = dcmplx(r1, r2)
        if (j == 0) Aref(il, iy, j) = Aref(il, iy, j) + 4.0d0      ! diagonally dominant
        if (iy + j >= ny) Aref(il, iy, j) = Aref(il, iy, j)*ph
        if (iy + j < 0) Aref(il, iy, j) = Aref(il, iy, j)*conjg(ph)
      end do
      xtrue(il, iy) = dcmplx(uniform_from_key(3, 1, iy, il, 0) - 0.5d0, uniform_from_key(4, 1, iy, il, 0) - 0.5d0)
    end do
  end do
  ! b = A x with columns taken modulo ny
  do il = 1, nl
    do iy = 0, ny - 1
      b = 0.0d0
      do j = -2, 2
        c = modulo(iy + j, ny)
        b = b + Aref(il, iy, j)*xtrue(il, c)
      end do
      X(il, iy) = b
      A(il, iy, :) = Aref(il, iy, :)
    end do
  end do
  !$omp target update to(A, X)
  call solve_lines(int(nl, C_INT))
  !$omp target update from(X)

  err = maxval(abs(X(1:nl, :) - xtrue))
  xmax = maxval(abs(xtrue))
  if (has_terminal) write (*, '(A,I0,A,I0,A,ES10.2,A,ES10.2)') '   cyclic pentadiagonal solve, ', nl, &
    ' lines of ', ny, ': max error ', err, '  (max |x| = ', xmax, ')'
  call free_linsolve()
  call free_mpi()
  call MPI_Finalize(ierr)
  if (err > 1.0d-12*xmax) then
    print *, 'FAILED'
    error stop 1
  end if
  if (has_terminal) print *, 'PASSED'

end program test_linsolve
