! The wall-normal (y) linear systems: one cyclic pentadiagonal system per
! Fourier mode, all modes solved at once, one GPU thread per mode.
!
! The five-point compact stencil on a periodic y gives a pentadiagonal
! matrix whose stencil wraps around the box; the wrapped entries carry the
! shear-periodic phase.  Written with band index j = -2..2 (column = row + j
! modulo ny), the matrix of a line is A(iy, j).
!
! The solve uses bordering: the last two unknowns x(ny-2), x(ny-1) are the
! border.  Rows and columns 0..ny-3 form a *plain* pentadiagonal matrix P
! (the wrapped entries of rows 0, 1 point at the border columns, and rows
! ny-2, ny-1 are the border rows), so
!   1. P is factored once, in place,
!   2. three right-hand sides are back-substituted: b, and the two columns
!      that couple to the border,
!   3. a 2x2 Schur complement gives the border, and the interior is corrected.
! This is the npy = 1 case of a distributed Schur solve.
!
! Storage is interleaved, line index first (A(iline, iy, j)), so that the
! threads of a kernel read consecutive addresses.  Lines are processed in
! batches of line_chunk x columns to bound the workspace.
module hst_linsolve

  use, intrinsic :: iso_c_binding
  use hst_params

  implicit none
  private
  public :: init_linsolve, free_linsolve, solve_lines, solve_component, apply_dy
  public :: A, X, Y1, Y2, nlines_max
  public :: KIND_D2V, KIND_ETA, KIND_D0

  integer(C_INT), parameter :: line_chunk = 16
  integer(C_INT), parameter :: KIND_D2V = 1, KIND_ETA = 2, KIND_D0 = 3
  complex(C_DOUBLE_COMPLEX), allocatable, save :: A(:, :, :), X(:, :), Y1(:, :), Y2(:, :)
  integer(C_INT), save :: nlines_max

contains

  subroutine init_linsolve()
    nlines_max = (2*nz + 1)*min(line_chunk, nxB)
    allocate (A(nlines_max, 0:ny - 1, -2:2), X(nlines_max, 0:ny - 1))
    allocate (Y1(nlines_max, 0:ny - 1), Y2(nlines_max, 0:ny - 1))
    A = 0; X = 0; Y1 = 0; Y2 = 0
    !$omp target enter data map(to: A, X, Y1, Y2)
  end subroutine init_linsolve

  subroutine free_linsolve()
    !$omp target exit data map(delete: A, X, Y1, Y2)
    deallocate (A, X, Y1, Y2)
  end subroutine free_linsolve

  ! Solve the systems of lines 1..nl held in A (matrices) and X (right-hand
  ! sides, overwritten by the solutions).  A is destroyed.
  subroutine solve_lines(nl)
    integer(C_INT), intent(in) :: nl
    integer(C_INT) :: il, n, ld
    n = ny
    ld = nlines_max
    !$omp target teams distribute parallel do default(none) shared(A, X, Y1, Y2, n, nl, ld) private(il)
    do il = 1, nl
      call cyclic_penta_solve(n, ld, il, A, X, Y1, Y2)
    end do
  end subroutine solve_lines

  ! One line.  Must stay inside this module: a declare-target procedure called
  ! across a module boundary does not survive nvlink.
  subroutine cyclic_penta_solve(n, ld, il, A, X, Y1, Y2)
    !$omp declare target
    integer(C_INT), intent(in) :: n, ld, il
    complex(C_DOUBLE_COMPLEX), intent(inout) :: A(ld, 0:n - 1, -2:2), X(ld, 0:n - 1)
    complex(C_DOUBLE_COMPLEX), intent(inout) :: Y1(ld, 0:n - 1), Y2(ld, 0:n - 1)
    integer(C_INT) :: i, m
    complex(C_DOUBLE_COMPLEX) :: l1, l2, d
    complex(C_DOUBLE_COMPLEX) :: c1m4, c1m3, c10, c2m3, c20, c21, d11, d12, d21, d22
    complex(C_DOUBLE_COMPLEX) :: s1, s2, m11, m12, m21, m22, det, xb1, xb2

    m = n - 2
    ! Columns n-2 (Y1) and n-1 (Y2) of the interior rows, then removed from P.
    do i = 0, m - 1
      Y1(il, i) = 0.0d0
      Y2(il, i) = 0.0d0
    end do
    Y1(il, 0) = A(il, 0, -2);      A(il, 0, -2) = 0.0d0
    Y2(il, 0) = A(il, 0, -1);      A(il, 0, -1) = 0.0d0
    Y2(il, 1) = A(il, 1, -2);      A(il, 1, -2) = 0.0d0
    Y1(il, m - 2) = A(il, m - 2, 2); A(il, m - 2, 2) = 0.0d0
    Y1(il, m - 1) = A(il, m - 1, 1); A(il, m - 1, 1) = 0.0d0
    Y2(il, m - 1) = A(il, m - 1, 2); A(il, m - 1, 2) = 0.0d0
    ! Border rows: their interior couplings (C) and their 2x2 block (D).
    c1m4 = A(il, n - 2, -2); c1m3 = A(il, n - 2, -1); c10 = A(il, n - 2, 2)
    c2m3 = A(il, n - 1, -2); c20 = A(il, n - 1, 1);   c21 = A(il, n - 1, 2)
    d11 = A(il, n - 2, 0); d12 = A(il, n - 2, 1)
    d21 = A(il, n - 1, -1); d22 = A(il, n - 1, 0)

    ! Factor P in place: multipliers overwrite the sub-diagonals.
    do i = 0, m - 2
      d = A(il, i, 0)
      l1 = A(il, i + 1, -1)/d
      A(il, i + 1, -1) = l1
      A(il, i + 1, 0) = A(il, i + 1, 0) - l1*A(il, i, 1)
      A(il, i + 1, 1) = A(il, i + 1, 1) - l1*A(il, i, 2)
      if (i + 2 <= m - 1) then
        l2 = A(il, i + 2, -2)/d
        A(il, i + 2, -2) = l2
        A(il, i + 2, -1) = A(il, i + 2, -1) - l2*A(il, i, 1)
        A(il, i + 2, 0) = A(il, i + 2, 0) - l2*A(il, i, 2)
      end if
    end do
    ! Three right-hand sides through the same factorisation.
    call penta_substitute(n, ld, il, m, A, X)
    call penta_substitute(n, ld, il, m, A, Y1)
    call penta_substitute(n, ld, il, m, A, Y2)
    ! 2x2 Schur complement for the border.
    s1 = X(il, n - 2) - (c1m4*X(il, m - 2) + c1m3*X(il, m - 1) + c10*X(il, 0))
    s2 = X(il, n - 1) - (c2m3*X(il, m - 1) + c20*X(il, 0) + c21*X(il, 1))
    m11 = d11 - (c1m4*Y1(il, m - 2) + c1m3*Y1(il, m - 1) + c10*Y1(il, 0))
    m12 = d12 - (c1m4*Y2(il, m - 2) + c1m3*Y2(il, m - 1) + c10*Y2(il, 0))
    m21 = d21 - (c2m3*Y1(il, m - 1) + c20*Y1(il, 0) + c21*Y1(il, 1))
    m22 = d22 - (c2m3*Y2(il, m - 1) + c20*Y2(il, 0) + c21*Y2(il, 1))
    det = m11*m22 - m12*m21
    xb1 = (m22*s1 - m12*s2)/det
    xb2 = (m11*s2 - m21*s1)/det
    X(il, n - 2) = xb1
    X(il, n - 1) = xb2
    do i = 0, m - 1
      X(il, i) = X(il, i) - Y1(il, i)*xb1 - Y2(il, i)*xb2
    end do
  end subroutine cyclic_penta_solve

  ! Forward and backward substitution through the factored interior P
  ! (rows 0..m-1) for one right-hand side b.
  subroutine penta_substitute(n, ld, il, m, A, b)
    !$omp declare target
    integer(C_INT), intent(in) :: n, ld, il, m
    complex(C_DOUBLE_COMPLEX), intent(in) :: A(ld, 0:n - 1, -2:2)
    complex(C_DOUBLE_COMPLEX), intent(inout) :: b(ld, 0:n - 1)
    integer(C_INT) :: k
    do k = 0, m - 2
      b(il, k + 1) = b(il, k + 1) - A(il, k + 1, -1)*b(il, k)
      if (k + 2 <= m - 1) b(il, k + 2) = b(il, k + 2) - A(il, k + 2, -2)*b(il, k)
    end do
    b(il, m - 1) = b(il, m - 1)/A(il, m - 1, 0)
    b(il, m - 2) = (b(il, m - 2) - A(il, m - 2, 1)*b(il, m - 1))/A(il, m - 2, 0)
    do k = m - 3, 0, -1
      b(il, k) = (b(il, k) - A(il, k, 1)*b(il, k + 1) - A(il, k, 2)*b(il, k + 2))/A(il, k, 0)
    end do
  end subroutine penta_substitute

  ! Assemble and solve one implicit system for component comp of V, whose
  ! interior rows hold the right-hand side on entry and the solution on exit:
  !   KIND_D2V :  lambda (D2 - k2 D0) - ni (D4 - 2 k2 D2 + k2^2 D0)
  !   KIND_ETA :  lambda D0 - ni (D2 - k2 D0)
  ! The (0,0) mode of KIND_D2V is singular and is set to zero.
  subroutine solve_component(kind, lambda, comp)
    integer(C_INT), intent(in) :: kind, comp
    real(C_DOUBLE), intent(in) :: lambda
    integer(C_INT) :: ix0, ix1, nl, ix, iz, iy, j, il, ncol
    real(C_DOUBLE) :: shift, coef
    complex(C_DOUBLE_COMPLEX) :: ph, wrap

    shift = modulo(S*time*ly, lx)
    do ix0 = nx0, nxN, line_chunk
      ix1 = min(ix0 + line_chunk - 1, nxN)
      ncol = 2*nz + 1
      nl = (ix1 - ix0 + 1)*ncol
      !$omp target teams distribute parallel do collapse(3) default(none) &
      !$omp shared(A, X, V, der, k2, ni, lambda, kind, comp, ix0, ix1, nz, ny, ncol, alfa0, shift) &
      !$omp private(ix, iz, iy, j, il, ph, wrap, coef)
      do ix = ix0, ix1
        do iz = -nz, nz
          do iy = 0, ny - 1
            il = (iz + nz + 1) + ncol*(ix - ix0)
            ph = exp(dcmplx(0.0d0, -alfa0*ix*shift))
            do j = -2, 2
              if (kind == KIND_D2V) then
                coef = lambda*(der(iy, 2, j) - k2(iz, ix)*der(iy, 0, j)) - &
                       ni*(der(iy, 3, j) - 2.0d0*k2(iz, ix)*der(iy, 2, j) + k2(iz, ix)*k2(iz, ix)*der(iy, 0, j))
              else
                coef = lambda*der(iy, 0, j) - ni*(der(iy, 2, j) - k2(iz, ix)*der(iy, 0, j))
              end if
              wrap = 1.0d0
              if (iy + j >= ny) wrap = ph
              if (iy + j < 0) wrap = conjg(ph)
              A(il, iy, j) = coef*wrap
            end do
            X(il, iy) = V(iy, iz, ix, comp)
          end do
        end do
      end do
      call solve_lines(nl)
      !$omp target teams distribute parallel do collapse(3) default(none) &
      !$omp shared(X, V, kind, comp, ix0, ix1, nz, ny, ncol) private(ix, iz, iy, il)
      do ix = ix0, ix1
        do iz = -nz, nz
          do iy = 0, ny - 1
            il = (iz + nz + 1) + ncol*(ix - ix0)
            if (kind == KIND_D2V .and. ix == 0 .and. iz == 0) then
              V(iy, iz, ix, comp) = 0.0d0
            else
              V(iy, iz, ix, comp) = X(il, iy)
            end if
          end do
        end do
      end do
    end do
  end subroutine solve_component

  ! dst = d/dy of V(:, :, :, src):  D0 x = D1 f, with the ghost rows of src
  ! already filled.  dst has the layout of one component of V (pass e.g.
  ! V(:, :, :, 3) or memrhs(:, :, :, 1)); only its interior rows are written.
  subroutine apply_dy(src, dst)
    integer(C_INT), intent(in) :: src
    complex(C_DOUBLE_COMPLEX), intent(inout) :: dst(ny0 - 2:, -nz:, nx0:)
    integer(C_INT) :: ix0, ix1, nl, ix, iz, iy, j, il, ncol
    real(C_DOUBLE) :: shift
    complex(C_DOUBLE_COMPLEX) :: ph, wrap, rhs

    shift = modulo(S*time*ly, lx)
    do ix0 = nx0, nxN, line_chunk
      ix1 = min(ix0 + line_chunk - 1, nxN)
      ncol = 2*nz + 1
      nl = (ix1 - ix0 + 1)*ncol
      !$omp target teams distribute parallel do collapse(3) default(none) &
      !$omp shared(A, X, V, der, src, ix0, ix1, nz, ny, ncol, alfa0, shift) &
      !$omp private(ix, iz, iy, j, il, ph, wrap, rhs)
      do ix = ix0, ix1
        do iz = -nz, nz
          do iy = 0, ny - 1
            il = (iz + nz + 1) + ncol*(ix - ix0)
            ph = exp(dcmplx(0.0d0, -alfa0*ix*shift))
            rhs = 0.0d0
            do j = -2, 2
              wrap = 1.0d0
              if (iy + j >= ny) wrap = ph
              if (iy + j < 0) wrap = conjg(ph)
              A(il, iy, j) = der(iy, 0, j)*wrap
              rhs = rhs + der(iy, 1, j)*V(iy + j, iz, ix, src)
            end do
            X(il, iy) = rhs
          end do
        end do
      end do
      call solve_lines(nl)
      !$omp target teams distribute parallel do collapse(3) default(none) &
      !$omp shared(X, dst, ix0, ix1, nz, ny, ncol) private(ix, iz, iy, il)
      do ix = ix0, ix1
        do iz = -nz, nz
          do iy = 0, ny - 1
            il = (iz + nz + 1) + ncol*(ix - ix0)
            dst(iy, iz, ix) = X(il, iy)
          end do
        end do
      end do
    end do
  end subroutine apply_dy

end module hst_linsolve
