! The wall-normal (y) linear systems: one cyclic pentadiagonal system per
! Fourier mode, all modes solved at once, one GPU thread per mode.
!
! The five-point compact stencil on a periodic y gives a pentadiagonal
! matrix whose stencil wraps around the box; the wrapped entries carry the
! shear-periodic phase.  Written with band index j = -2..2 (column = row + j
! modulo ny), row iy of a line is
!   a(j) = coef(kind, iy, j) * wrap,   wrap = ph if iy + j >= ny, conjg(ph) if iy + j < 0, else 1
! with the real stencil combination coef of the system kind (line_solve).
!
! The solve uses bordering: the last two unknowns x(ny-2), x(ny-1) are the
! border.  Rows and columns 0..ny-3 form a *plain* pentadiagonal matrix P
! (the wrapped entries of rows 0, 1 point at the border columns, and rows
! ny-2, ny-1 are the border rows), so
!   1. P is factored, and three right-hand sides are forward-substituted in
!      the same sweep: b, and the two columns that couple to the border,
!   2. the three are back-substituted,
!   3. a 2x2 Schur complement gives the border, and the interior is corrected.
! This is the npy = 1 case of a distributed Schur solve.
!
! The rows are generated on the fly inside the sweep, so the matrix is
! never stored: the forward sweep keeps the two previous rows in registers
! and writes only the three upper diagonals U (for the back substitution)
! and the substituted right-hand sides X, Y1, Y2.  Storage is interleaved,
! line index first (U(iline, iy, j)), so that the threads of a kernel read
! consecutive addresses.  Lines are processed in batches of line_chunk x
! columns (deck parameter; 0 = all columns at once on the GPU, 16 on the
! CPU) to bound that workspace.
module hst_linsolve

  use, intrinsic :: iso_c_binding
  use hst_params
  use hst_derivatives, only: shear_shifts

  implicit none
  private
  public :: init_linsolve, free_linsolve, line_solve
  public :: KIND_D2V, KIND_ETA, KIND_POISSON, KIND_D0, KIND_DY

  ! the systems line_solve assembles (see there)
  integer(C_INT), parameter :: KIND_D2V = 1, KIND_ETA = 2, KIND_POISSON = 3, KIND_D0 = 4, KIND_DY = 5
  complex(C_DOUBLE_COMPLEX), allocatable, save :: U(:, :, :), X(:, :), Y1(:, :), Y2(:, :)
  integer(C_INT), save :: nlines_max, chunk

contains

  subroutine init_linsolve()
    ! default: every column at once on the GPU (fewest launches); 16 on the
    ! CPU, where the workspace of all columns falls out of the cache
    chunk = 16
#ifdef HAVE_CUDA
    chunk = nxB
#endif
    if (line_chunk > 0) chunk = line_chunk
    chunk = min(chunk, nxB)
    nlines_max = (2*nz + 1)*chunk
    allocate (U(nlines_max, 0:ny - 1, 0:2), X(nlines_max, 0:ny - 1))
    allocate (Y1(nlines_max, 0:ny - 1), Y2(nlines_max, 0:ny - 1))
    U = 0; X = 0; Y1 = 0; Y2 = 0
    !$omp target enter data map(to: U, X, Y1, Y2)
    if (has_terminal) write (*, '(A,I0,A,I0,A,F8.1,A)') '   line solver: batches of ', chunk, ' x columns, ', &
      nlines_max, ' lines, workspace ', 6.0d0*16.0d0*nlines_max*ny/1024.0d0**2, ' MB per rank'
  end subroutine init_linsolve

  subroutine free_linsolve()
    !$omp target exit data map(delete: U, X, Y1, Y2)
    deallocate (U, X, Y1, Y2)
  end subroutine free_linsolve

  ! Assemble and solve one system per mode for fields with the layout of a
  ! component of V (pass e.g. V(:, :, :, 2) or rhs(:, :, :, 1)); the
  ! right-hand side comes from src and the solution goes to the interior
  ! rows of dst (a different array):
  !   KIND_D2V     [lambda (D2 - k2 D0) - ni (D4 - 2 k2 D2 + k2^2 D0)] x = src
  !   KIND_ETA     [lambda D0 - ni (D2 - k2 D0)] x = src
  !   KIND_POISSON (D2 - k2 D0) x = src                  (lambda unused)
  !   KIND_D0      D0 x = src         the unweighted quantity behind a D0-weighted sum
  !   KIND_DY      D0 x = D1 src      d/dy, with the ghost rows of src filled
  ! The wrap phase uses the displacements shift_x, shift_z of the upper
  ! image; without them, those of the current time.  The (0,0) mode of
  ! KIND_D2V is singular and is set to zero; that of KIND_POISSON is
  ! singular too and dst is left untouched there for the caller.
  subroutine line_solve(kind, lambda, src, dst, shift_x, shift_z)
    integer(C_INT), intent(in) :: kind
    real(C_DOUBLE), intent(in) :: lambda
    complex(C_DOUBLE_COMPLEX), intent(in) :: src(ny0 - 2:, -nz:, nx0:)
    complex(C_DOUBLE_COMPLEX), intent(inout) :: dst(ny0 - 2:, -nz:, nx0:)
    real(C_DOUBLE), intent(in), optional :: shift_x, shift_z
    integer(C_INT) :: ix0, ix1, nl, ix, iz, iy, j, il, ncol, n, ld
    real(C_DOUBLE) :: sx, sz
    complex(C_DOUBLE_COMPLEX) :: b, ph

    if (present(shift_x)) then
      sx = shift_x; sz = shift_z
    else
      call shear_shifts(time, sx, sz)
    end if
    ncol = 2*nz + 1
    n = ny
    ld = nlines_max
    do ix0 = nx0, nxN, chunk
      ix1 = min(ix0 + chunk - 1, nxN)
      nl = (ix1 - ix0 + 1)*ncol
      ! gather the right-hand sides into the interleaved layout
      !$omp target teams distribute parallel do collapse(3) default(none) &
      !$omp shared(X, src, der, kind, ix0, ix1, nz, ny, ncol) private(ix, iz, iy, j, il, b)
      do ix = ix0, ix1
        do iz = -nz, nz
          do iy = 0, ny - 1
            il = (iz + nz + 1) + ncol*(ix - ix0)
            if (kind == KIND_DY) then
              b = 0.0d0
              do j = -2, 2
                b = b + der(iy, 1, j)*src(iy + j, iz, ix)
              end do
            else
              b = src(iy, iz, ix)
            end if
            X(il, iy) = b
          end do
        end do
      end do
      ! one thread per line: build the rows, factor and solve
      !$omp target teams distribute parallel do default(none) &
      !$omp shared(U, X, Y1, Y2, der, k2, ni, lambda, kind, ix0, nz, ncol, alfa0, beta0, sx, sz, n, nl, ld) &
      !$omp private(il, ix, iz, ph)
      do il = 1, nl
        ix = ix0 + (il - 1)/ncol
        iz = mod(il - 1, ncol) - nz
        ph = exp(dcmplx(0.0d0, -(alfa0*ix*sx + beta0*iz*sz)))
        call cyclic_penta_solve(kind, lambda, ni, k2(iz, ix), ph, n, ld, il, der, U, X, Y1, Y2)
      end do
      ! scatter the solutions
      !$omp target teams distribute parallel do collapse(3) default(none) &
      !$omp shared(X, dst, kind, ix0, ix1, nz, ny, ncol) private(ix, iz, iy, il)
      do ix = ix0, ix1
        do iz = -nz, nz
          do iy = 0, ny - 1
            il = (iz + nz + 1) + ncol*(ix - ix0)
            if (ix == 0 .and. iz == 0 .and. kind == KIND_D2V) then
              dst(iy, iz, ix) = 0.0d0
            else if (.not. (ix == 0 .and. iz == 0 .and. kind == KIND_POISSON)) then
              dst(iy, iz, ix) = X(il, iy)
            end if
          end do
        end do
      end do
    end do
  end subroutine line_solve

  ! Row iy, band entry j of the system `kind` before the wrap phase.
  real(C_DOUBLE) function coef(kind, lambda, ni, kk, n, der, iy, j)
    !$omp declare target
    integer(C_INT), intent(in) :: kind, n, iy, j
    real(C_DOUBLE), intent(in) :: lambda, ni, kk
    real(C_DOUBLE), intent(in) :: der(0:n - 1, 0:3, -2:2)
    select case (kind)
    case (KIND_D2V)
      coef = lambda*(der(iy, 2, j) - kk*der(iy, 0, j)) - &
             ni*(der(iy, 3, j) - 2.0d0*kk*der(iy, 2, j) + kk*kk*der(iy, 0, j))
    case (KIND_ETA)
      coef = lambda*der(iy, 0, j) - ni*(der(iy, 2, j) - kk*der(iy, 0, j))
    case (KIND_POISSON)
      coef = der(iy, 2, j) - kk*der(iy, 0, j)
    case default
      coef = der(iy, 0, j)
    end select
  end function coef

  ! One line, right-hand side and solution in X(il, :).  Must stay inside
  ! this module: a declare-target procedure called across a module boundary
  ! does not survive nvlink.
  subroutine cyclic_penta_solve(kind, lambda, ni, kk, ph, n, ld, il, der, U, X, Y1, Y2)
    !$omp declare target
    integer(C_INT), intent(in) :: kind, n, ld, il
    real(C_DOUBLE), intent(in) :: lambda, ni, kk
    complex(C_DOUBLE_COMPLEX), intent(in) :: ph
    real(C_DOUBLE), intent(in) :: der(0:n - 1, 0:3, -2:2)
    complex(C_DOUBLE_COMPLEX), intent(inout) :: U(ld, 0:n - 1, 0:2), X(ld, 0:n - 1)
    complex(C_DOUBLE_COMPLEX), intent(inout) :: Y1(ld, 0:n - 1), Y2(ld, 0:n - 1)
    integer(C_INT) :: i, j, m
    complex(C_DOUBLE_COMPLEX) :: a(-2:2), wrap, l1, l2, bx, b1, b2          ! row i: substituted right-hand sides
    complex(C_DOUBLE_COMPLEX) :: p1, u11, u12, bx1, b11, b21      ! row i-1: pivot, upper entries, substituted rhs
    complex(C_DOUBLE_COMPLEX) :: p2, u21, u22, bx2, b12, b22      ! row i-2
    complex(C_DOUBLE_COMPLEX) :: c1m4, c1m3, c10, c2m3, c20, c21, d11, d12, d21, d22
    complex(C_DOUBLE_COMPLEX) :: s1, s2, m11, m12, m21, m22, det, xb1, xb2
    complex(C_DOUBLE_COMPLEX) :: xk1, xk2, yk1, yk2, zk1, zk2

    m = n - 2
    ! Forward sweep over the interior rows: generate row i, move its border
    ! columns (n-2 -> Y1, n-1 -> Y2) out of P, eliminate with rows i-2 and
    ! i-1, store the upper part and the substituted right-hand sides.
    p1 = 1.0d0; u11 = 0.0d0; u12 = 0.0d0; bx1 = 0.0d0; b11 = 0.0d0; b21 = 0.0d0
    p2 = 1.0d0; u21 = 0.0d0; u22 = 0.0d0; bx2 = 0.0d0; b12 = 0.0d0; b22 = 0.0d0
    do i = 0, m - 1
      do j = -2, 2
        wrap = 1.0d0
        if (i + j >= n) wrap = ph
        if (i + j < 0) wrap = conjg(ph)
        a(j) = coef(kind, lambda, ni, kk, n, der, i, j)*wrap
      end do
      b1 = 0.0d0; b2 = 0.0d0
      if (i == 0) then
        b1 = a(-2); b2 = a(-1); a(-2) = 0.0d0; a(-1) = 0.0d0
      end if
      if (i == 1) then
        b2 = a(-2); a(-2) = 0.0d0
      end if
      if (i == m - 2) then
        b1 = a(2); a(2) = 0.0d0
      end if
      if (i == m - 1) then
        b1 = a(1); b2 = a(2); a(1) = 0.0d0; a(2) = 0.0d0
      end if
      bx = X(il, i)
      if (i >= 2) then
        l2 = a(-2)/p2
        a(-1) = a(-1) - l2*u21
        a(0) = a(0) - l2*u22
        bx = bx - l2*bx2; b1 = b1 - l2*b12; b2 = b2 - l2*b22
      end if
      if (i >= 1) then
        l1 = a(-1)/p1
        a(0) = a(0) - l1*u11
        a(1) = a(1) - l1*u12
        bx = bx - l1*bx1; b1 = b1 - l1*b11; b2 = b2 - l1*b21
      end if
      U(il, i, 0) = a(0); U(il, i, 1) = a(1); U(il, i, 2) = a(2)
      X(il, i) = bx; Y1(il, i) = b1; Y2(il, i) = b2
      p2 = p1; u21 = u11; u22 = u12; bx2 = bx1; b12 = b11; b22 = b21
      p1 = a(0); u11 = a(1); u12 = a(2); bx1 = bx; b11 = b1; b21 = b2
    end do
    ! Border rows: their interior couplings (C) and their 2x2 block (D).
    c1m4 = coef(kind, lambda, ni, kk, n, der, n - 2, -2)
    c1m3 = coef(kind, lambda, ni, kk, n, der, n - 2, -1)
    c10 = coef(kind, lambda, ni, kk, n, der, n - 2, 2)*ph
    c2m3 = coef(kind, lambda, ni, kk, n, der, n - 1, -2)
    c20 = coef(kind, lambda, ni, kk, n, der, n - 1, 1)*ph
    c21 = coef(kind, lambda, ni, kk, n, der, n - 1, 2)*ph
    d11 = coef(kind, lambda, ni, kk, n, der, n - 2, 0)
    d12 = coef(kind, lambda, ni, kk, n, der, n - 2, 1)
    d21 = coef(kind, lambda, ni, kk, n, der, n - 1, -1)
    d22 = coef(kind, lambda, ni, kk, n, der, n - 1, 0)
    ! Backward sweep for the three right-hand sides.
    xk1 = X(il, m - 1)/U(il, m - 1, 0); yk1 = Y1(il, m - 1)/U(il, m - 1, 0); zk1 = Y2(il, m - 1)/U(il, m - 1, 0)
    X(il, m - 1) = xk1; Y1(il, m - 1) = yk1; Y2(il, m - 1) = zk1
    xk2 = xk1; yk2 = yk1; zk2 = zk1
    xk1 = (X(il, m - 2) - U(il, m - 2, 1)*xk1)/U(il, m - 2, 0)
    yk1 = (Y1(il, m - 2) - U(il, m - 2, 1)*yk1)/U(il, m - 2, 0)
    zk1 = (Y2(il, m - 2) - U(il, m - 2, 1)*zk1)/U(il, m - 2, 0)
    X(il, m - 2) = xk1; Y1(il, m - 2) = yk1; Y2(il, m - 2) = zk1
    do i = m - 3, 0, -1
      bx = (X(il, i) - U(il, i, 1)*xk1 - U(il, i, 2)*xk2)/U(il, i, 0)
      b1 = (Y1(il, i) - U(il, i, 1)*yk1 - U(il, i, 2)*yk2)/U(il, i, 0)
      b2 = (Y2(il, i) - U(il, i, 1)*zk1 - U(il, i, 2)*zk2)/U(il, i, 0)
      X(il, i) = bx; Y1(il, i) = b1; Y2(il, i) = b2
      xk2 = xk1; yk2 = yk1; zk2 = zk1
      xk1 = bx; yk1 = b1; zk1 = b2
    end do
    ! 2x2 Schur complement for the border, then the interior correction.
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

end module hst_linsolve
