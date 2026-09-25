! Unit test of the cyclic pentadiagonal line solver (WP2).
!
!   mpirun -np N build-cpu/test_linsolve [hst.in]
!
! For each system kind of line_solve the operator is applied on the host to
! a random field: the stencil sums run over the ghost rows, which hold the
! shear-periodic images at a nonzero time, so the wrap phase of the solver
! is checked against the ghost-row convention of the rest of the code.  The
! solve of that right-hand side must give the field back to round-off.
! KIND_DY is checked through its defining identity D0 (dst) = D1 src.
! Runs on the device in a GPU build.
program test_linsolve

  use, intrinsic :: iso_c_binding
  use mpi_f08
  use hst_params
  use test_common
  use hst_linsolve
  use hst_derivatives, only: shear_shifts, fill_ghosts
  use hst_initial, only: uniform_from_key

  implicit none

  integer :: ierr, ix, iy, iz, j, kind
  real(C_DOUBLE) :: lambda, kk, c, err, ref, worst, err_g, ref_g, sx, sz
  complex(C_DOUBLE_COMPLEX) :: acc, ph
  character(len=8), parameter :: name(5) = ['D2V     ', 'ETA     ', 'POISSON ', 'D0      ', 'DY      ']

  call test_start()
  time = 0.37d0                              ! a nonzero wrap phase
  call test_setup()
  call shear_shifts(time, sx, sz)

  ! random field in V(:, :, :, 1), (0,0) mode zero (singular for two kinds),
  ! ghost rows at the current time
  V = 0
  do ix = nx0, nxN
    do iz = -nz, nz
      if (ix == 0 .and. iz == 0) cycle
      do iy = 0, ny - 1
        V(iy, iz, ix, 1) = dcmplx(uniform_from_key(3, 1, iy, iz, ix) - 0.5d0, uniform_from_key(4, 1, iy, iz, ix) - 0.5d0)
      end do
    end do
  end do
  !$omp target update to(V)
  call fill_ghosts(1)
  !$omp target update from(V)
  lambda = 123.4d0
  worst = 0.0d0

  do kind = KIND_D2V, KIND_DY
    ! right-hand side on the host: the operator applied through the ghost rows
    do ix = nx0, nxN
      do iz = -nz, nz
        kk = k2(iz, ix)
        do iy = 0, ny - 1
          acc = 0.0d0
          do j = -2, 2
            select case (kind)
            case (KIND_D2V)
              c = lambda*(der(iy, 2, j) - kk*der(iy, 0, j)) - &
                  ni*(der(iy, 3, j) - 2.0d0*kk*der(iy, 2, j) + kk*kk*der(iy, 0, j))
            case (KIND_ETA)
              c = lambda*der(iy, 0, j) - ni*(der(iy, 2, j) - kk*der(iy, 0, j))
            case (KIND_POISSON)
              c = der(iy, 2, j) - kk*der(iy, 0, j)
            case default
              c = der(iy, 0, j)
            end select
            acc = acc + c*V(iy + j, iz, ix, 1)
          end do
          rhs(iy, iz, ix, 1) = acc
        end do
      end do
    end do
    rhs(:, :, :, 2) = 0
    !$omp target update to(rhs)
    if (kind == KIND_DY) then
      call line_solve(kind, lambda, V(:, :, :, 1), rhs(:, :, :, 2))
    else
      call line_solve(kind, lambda, rhs(:, :, :, 1), rhs(:, :, :, 2))
    end if
    !$omp target update from(rhs)
    err = 0.0d0; ref = 0.0d0
    if (kind == KIND_DY) then
      ! D0 dst must equal D1 src: fill the ghost rows of dst on the host
      do ix = nx0, nxN
        do iz = -nz, nz
          ph = exp(dcmplx(0.0d0, -(alfa0*ix*sx + beta0*iz*sz)))
          rhs(ny, iz, ix, 2) = rhs(0, iz, ix, 2)*ph; rhs(ny + 1, iz, ix, 2) = rhs(1, iz, ix, 2)*ph
          rhs(-1, iz, ix, 2) = rhs(ny - 1, iz, ix, 2)*conjg(ph); rhs(-2, iz, ix, 2) = rhs(ny - 2, iz, ix, 2)*conjg(ph)
          do iy = 0, ny - 1
            acc = 0.0d0
            do j = -2, 2
              acc = acc + der(iy, 0, j)*rhs(iy + j, iz, ix, 2) - der(iy, 1, j)*V(iy + j, iz, ix, 1)
            end do
            err = max(err, abs(acc))
            ref = max(ref, abs(rhs(iy, iz, ix, 1)))     ! the D0-weighted field itself
          end do
        end do
      end do
    else
      err = maxval(abs(rhs(0:ny - 1, :, :, 2) - V(0:ny - 1, :, :, 1)))
      ref = maxval(abs(V(0:ny - 1, :, :, 1)))
    end if
    call MPI_Allreduce(err, err_g, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(ref, ref_g, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    worst = max(worst, err_g/ref_g)
    if (has_terminal) write (*, '(A,A,A,ES10.2,A,ES10.2)') '   ', name(kind), ': max error ', err_g, '  relative ', err_g/ref_g
  end do
  if (has_terminal) write (*, '(A,I0,A,I0,A)') '   (', (2*nz + 1)*nxB, ' lines of ', ny, ' on rank 0)'
  call test_finish(worst <= 1.0d-11)

end program test_linsolve
