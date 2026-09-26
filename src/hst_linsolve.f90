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
!   2. a 2x2 Schur complement gives the border,
!   3. the back substitution of b minus the border columns times the border
!      is the solution.
! This is the npy = 1 case of a distributed Schur solve.
!
! The rows are generated on the fly inside the sweep, so the matrix is
! never stored: the forward sweep keeps the two previous rows in registers,
! divides each row by its pivot (one complex division per row) and writes
! only the two upper entries U1, U2 (for the back substitution) and the
! substituted right-hand sides X, Y1, Y2.  The Schur complement needs the
! back-substituted values at the four rows coupled to the border: the last
! two are the last two rows of the sweep, and the first two are inner
! products of the substituted right-hand sides with the first two rows of
! the inverse of the unit upper factor, which is a *forward* recurrence
! (w(i) = -U1(i-1) w(i-1) - U2(i-2) w(i-2)), so both are accumulated in
! the same sweep and the border is known before the back substitution.
! The solve is then two passes over the line: the forward sweep reads the
! right-hand side from the field and writes the five columns, the backward
! sweep reads them and writes the solution into the field (192 bytes per
! row instead of the 384 of a gather, a three-pass solve and a scatter).
! Storage is interleaved, line index first (U1(iline, iy)), so that the
! threads of a kernel read consecutive addresses.  Lines are processed in
! batches of line_chunk x columns (deck parameter; 0 = all columns at once
! on the GPU, 16 on the CPU) to bound that workspace.
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
  complex(C_DOUBLE_COMPLEX), allocatable, save :: U1(:, :), U2(:, :), X(:, :), Y1(:, :), Y2(:, :)
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
    allocate (U1(nlines_max, 0:ny - 1), U2(nlines_max, 0:ny - 1), X(nlines_max, 0:ny - 1))
    allocate (Y1(nlines_max, 0:ny - 1), Y2(nlines_max, 0:ny - 1))
    U1 = 0; U2 = 0; X = 0; Y1 = 0; Y2 = 0
    !$omp target enter data map(to: U1, U2, X, Y1, Y2)
    if (has_terminal) write (*, '(A,I0,A,I0,A,F8.1,A)') '   line solver: batches of ', chunk, ' x columns, ', &
      nlines_max, ' lines, workspace ', 5.0d0*16.0d0*nlines_max*ny/1024.0d0**2, ' MB per rank'
  end subroutine init_linsolve

  subroutine free_linsolve()
    !$omp target exit data map(delete: U1, U2, X, Y1, Y2)
    deallocate (U1, U2, X, Y1, Y2)
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
    complex(C_DOUBLE_COMPLEX), intent(in), contiguous :: src(ny0 - 2:, -nz:, nx0:)
    complex(C_DOUBLE_COMPLEX), intent(inout), contiguous :: dst(ny0 - 2:, -nz:, nx0:)
    real(C_DOUBLE), intent(in), optional :: shift_x, shift_z
    integer(C_INT) :: ix0, ix1, nl, ix, iz, iy, il, ncol, n, ld
    real(C_DOUBLE) :: sx, sz
    complex(C_DOUBLE_COMPLEX) :: ph

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
      ! one thread per line: build the rows, factor and solve
      !$omp target teams distribute parallel do default(none) &
      !$omp shared(U1, U2, X, Y1, Y2, src, dst, der, k2, ni, lambda, kind, ix0, nz, nx0, ncol, alfa0, beta0, sx, sz, n, nl, ld) &
      !$omp private(il, ix, iz, iy, ph)
      do il = 1, nl
        ix = ix0 + (il - 1)/ncol
        iz = mod(il - 1, ncol) - nz
        if (ix == 0 .and. iz == 0 .and. kind == KIND_D2V) then
          do iy = 0, n - 1
            dst(iy, iz, ix) = 0.0d0
          end do
        else if (.not. (ix == 0 .and. iz == 0 .and. kind == KIND_POISSON)) then
          ph = exp(dcmplx(0.0d0, -(alfa0*ix*sx + beta0*iz*sz)))
          call cyclic_penta_solve(kind, lambda, ni, k2(iz, ix), ph, n, ld, il, (iz + nz + 1) + ncol*(ix - nx0), &
                                  der, U1, U2, X, Y1, Y2, src, dst)
        end if
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

  ! Right-hand side of row iy of line jl of src (the field's lines in memory
  ! order): the field itself, or its D1 stencil through the ghost rows for
  ! KIND_DY.
  complex(C_DOUBLE_COMPLEX) function rhs_row(kind, n, der, src, iy, jl)
    !$omp declare target
    integer(C_INT), intent(in) :: kind, n, iy, jl
    real(C_DOUBLE), intent(in) :: der(0:n - 1, 0:3, -2:2)
    complex(C_DOUBLE_COMPLEX), intent(in) :: src(-2:n + 1, *)
    integer(C_INT) :: j
    if (kind == KIND_DY) then
      rhs_row = 0.0d0
      do j = -2, 2
        rhs_row = rhs_row + der(iy, 1, j)*src(iy + j, jl)
      end do
    else
      rhs_row = src(iy, jl)
    end if
  end function rhs_row

  ! One line: right-hand side in src(:, jl), solution to dst(:, jl) (line jl
  ! of the field in memory order, ghost rows included), workspace column il.
  ! Must stay inside this module: a declare-target procedure called across
  ! a module boundary does not survive nvlink.
  subroutine cyclic_penta_solve(kind, lambda, ni, kk, ph, n, ld, il, jl, der, U1, U2, X, Y1, Y2, src, dst)
    !$omp declare target
    integer(C_INT), intent(in) :: kind, n, ld, il, jl
    real(C_DOUBLE), intent(in) :: lambda, ni, kk
    complex(C_DOUBLE_COMPLEX), intent(in) :: ph
    real(C_DOUBLE), intent(in) :: der(0:n - 1, 0:3, -2:2)
    complex(C_DOUBLE_COMPLEX), intent(inout) :: U1(ld, 0:n - 1), U2(ld, 0:n - 1), X(ld, 0:n - 1)
    complex(C_DOUBLE_COMPLEX), intent(inout) :: Y1(ld, 0:n - 1), Y2(ld, 0:n - 1)
    complex(C_DOUBLE_COMPLEX), intent(in) :: src(-2:n + 1, *)
    complex(C_DOUBLE_COMPLEX), intent(inout) :: dst(-2:n + 1, *)
    integer(C_INT) :: i, j, m
    complex(C_DOUBLE_COMPLEX) :: a(-2:2), wrap, l, rp, t1, t2, bx, p, q     ! row i: upper entries and right-hand sides (b, border columns), scaled by 1/pivot
    complex(C_DOUBLE_COMPLEX) :: u11, u12, x1, p1, q1                       ! row i-1, scaled
    complex(C_DOUBLE_COMPLEX) :: u21, u22, x2, p2, q2                       ! row i-2, scaled
    complex(C_DOUBLE_COMPLEX) :: w, w1, w2, v, v1, v2                       ! rows 0 and 1 of the inverse unit upper factor at i, i-1, i-2
    complex(C_DOUBLE_COMPLEX) :: sx0, sp0, sq0, sx1, sp1, sq1               ! their inner products with the substituted right-hand sides
    complex(C_DOUBLE_COMPLEX) :: xm1, pm1, qm1, xm2, pm2, qm2               ! back-substituted rows m-1, m-2
    complex(C_DOUBLE_COMPLEX) :: c1m4, c1m3, c10, c2m3, c20, c21, d11, d12, d21, d22
    complex(C_DOUBLE_COMPLEX) :: s1, s2, m11, m12, m21, m22, det, xb1, xb2, xk, xk1, xk2

    m = n - 2
    ! Forward sweep over the interior rows: generate row i, move its border
    ! columns (n-2 -> p, n-1 -> q) out of P, eliminate with rows i-2 and
    ! i-1, divide by the pivot, store the upper part and the substituted
    ! right-hand sides; accumulate the first two rows of the back substitution.
    u11 = 0.0d0; u12 = 0.0d0; x1 = 0.0d0; p1 = 0.0d0; q1 = 0.0d0
    u21 = 0.0d0; u22 = 0.0d0; x2 = 0.0d0; p2 = 0.0d0; q2 = 0.0d0
    w1 = 0.0d0; w2 = 0.0d0; v1 = 0.0d0; v2 = 0.0d0
    sx0 = 0.0d0; sp0 = 0.0d0; sq0 = 0.0d0; sx1 = 0.0d0; sp1 = 0.0d0; sq1 = 0.0d0
    do i = 0, m - 1
      do j = -2, 2
        wrap = 1.0d0
        if (i + j >= n) wrap = ph
        if (i + j < 0) wrap = conjg(ph)
        a(j) = coef(kind, lambda, ni, kk, n, der, i, j)*wrap
      end do
      p = 0.0d0; q = 0.0d0
      if (i == 0) then
        p = a(-2); q = a(-1); a(-2) = 0.0d0; a(-1) = 0.0d0
      end if
      if (i == 1) then
        q = a(-2); a(-2) = 0.0d0
      end if
      if (i == m - 2) then
        p = a(2); a(2) = 0.0d0
      end if
      if (i == m - 1) then
        p = a(1); q = a(2); a(1) = 0.0d0; a(2) = 0.0d0
      end if
      bx = rhs_row(kind, n, der, src, i, jl)
      if (i >= 2) then
        l = a(-2)
        a(-1) = a(-1) - l*u21
        a(0) = a(0) - l*u22
        bx = bx - l*x2; p = p - l*p2; q = q - l*q2
      end if
      if (i >= 1) then
        l = a(-1)
        a(0) = a(0) - l*u11
        a(1) = a(1) - l*u12
        bx = bx - l*x1; p = p - l*p1; q = q - l*q1
      end if
      rp = 1.0d0/a(0)
      t1 = a(1)*rp; t2 = a(2)*rp; bx = bx*rp; p = p*rp; q = q*rp
      U1(il, i) = t1; U2(il, i) = t2; X(il, i) = bx; Y1(il, i) = p; Y2(il, i) = q
      ! rows 0 and 1 of the inverse of the unit upper factor: w(0) = 1, v(1) = 1,
      ! w(i) = -U1(i-1) w(i-1) - U2(i-2) w(i-2), the same for v
      w = -(u11*w1 + u22*w2); if (i == 0) w = 1.0d0
      v = -(u11*v1 + u22*v2); if (i == 1) v = 1.0d0
      sx0 = sx0 + w*bx; sp0 = sp0 + w*p; sq0 = sq0 + w*q
      sx1 = sx1 + v*bx; sp1 = sp1 + v*p; sq1 = sq1 + v*q
      w2 = w1; w1 = w; v2 = v1; v1 = v
      u21 = u11; u22 = u12; x2 = x1; p2 = p1; q2 = q1
      u11 = t1; u12 = t2; x1 = bx; p1 = p; q1 = q
    end do
    ! Back-substituted values at the rows coupled to the border: m-1 (its
    ! upper entries are border columns), m-2 (U2 = 0 there), and 0, 1 (the sums).
    xm1 = x1; pm1 = p1; qm1 = q1
    xm2 = x2 - u21*xm1; pm2 = p2 - u21*pm1; qm2 = q2 - u21*qm1
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
    ! 2x2 Schur complement for the border.
    s1 = rhs_row(kind, n, der, src, n - 2, jl) - (c1m4*xm2 + c1m3*xm1 + c10*sx0)
    s2 = rhs_row(kind, n, der, src, n - 1, jl) - (c2m3*xm1 + c20*sx0 + c21*sx1)
    m11 = d11 - (c1m4*pm2 + c1m3*pm1 + c10*sp0)
    m12 = d12 - (c1m4*qm2 + c1m3*qm1 + c10*sq0)
    m21 = d21 - (c2m3*pm1 + c20*sp0 + c21*sp1)
    m22 = d22 - (c2m3*qm1 + c20*sq0 + c21*sq1)
    det = m11*m22 - m12*m21
    xb1 = (m22*s1 - m12*s2)/det
    xb2 = (m11*s2 - m21*s1)/det
    dst(n - 2, jl) = xb1
    dst(n - 1, jl) = xb2
    ! Backward sweep of b minus the border columns times the border, into dst
    ! (U1(m-1), U2(m-1), U2(m-2) are zero, so one recurrence serves all rows).
    xk1 = 0.0d0; xk2 = 0.0d0
    do i = m - 1, 0, -1
      xk = X(il, i) - Y1(il, i)*xb1 - Y2(il, i)*xb2 - U1(il, i)*xk1 - U2(il, i)*xk2
      dst(i, jl) = xk
      xk2 = xk1; xk1 = xk
    end do
  end subroutine cyclic_penta_solve

end module hst_linsolve
