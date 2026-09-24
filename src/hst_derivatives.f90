! Compact finite differences in y, and the shear-periodic ghost rows.
!
! der(iy, n, j), j = -2..2, are the five-point weights at node iy for
!   n = 0 : the compact "interpolation" operator D0,  1 : D1,  2 : D2,  3 : D4,
! chosen so that the scheme  sum_j der(iy,0,j) f'(iy+j) = sum_j der(iy,1,j) f(iy+j)
! (and likewise for D2, D4) is exact on polynomials.  Every equation in the
! code is written in this D0-weighted form, so a derivative is a solve with
! D0 on the left.  From channel/src/run/case_setup.fypp setup_derivatives,
! which is the CPL setup_derivatives; the one-sided wall rows are gone.
!
! Ghost rows: the box is periodic in y up to the mean-shear displacement,
!   f(x, y + ly) = f(x - S t ly, y)   ->   f_hat(ny) = f_hat(0) * exp(-i kx S t ly)
! so the two rows above the box are the two lowest rows times that phase and
! the two rows below are the two highest rows times its conjugate.  With the
! images in place every stencil in the code is a plain five-point sum.
module hst_derivatives

  use, intrinsic :: iso_c_binding
  use hst_params

  implicit none
  private
  public :: setup_derivatives, fill_ghosts, shear_shift_length

contains

  subroutine setup_derivatives()
    real(C_DOUBLE) :: M(0:4, 0:4), t(0:4), dyj(-2:2)
    integer(C_INT) :: iy, i, j

    if (ny < 8) error stop 'ny must be at least 8'
    allocate (der(0:ny - 1, 0:3, -2:2))
    do iy = 0, ny - 1
      dyj = y(iy - 2:iy + 2) - y(iy)
      ! D4: exact fourth derivative of y^4 (= 24) and zero on y^0..y^3
      do i = 0, 4; do j = 0, 4; M(i, j) = dyj(j - 2)**(4 - i); end do; end do
      t = 0; t(0) = 24
      der(iy, 3, :) = solve5(M, t)
      ! D0: sum_j d0_j f''''(y_j) = sum_j d4_j f(y_j) exact for f = y^8 .. y^4
      do i = 0, 4; do j = 0, 4
        M(i, j) = (5 - i)*(6 - i)*(7 - i)*(8 - i)*dyj(j - 2)**(4 - i)
      end do; end do
      do i = 0, 4; t(i) = sum(der(iy, 3, :)*dyj**(8 - i)); end do
      der(iy, 0, :) = solve5(M, t)
      ! D2 and D1 consistent with that D0
      do i = 0, 4; do j = 0, 4; M(i, j) = dyj(j - 2)**(4 - i); end do; end do
      t = 0; do i = 0, 2; t(i) = sum(der(iy, 0, :)*(4 - i)*(3 - i)*dyj**(2 - i)); end do
      der(iy, 2, :) = solve5(M, t)
      t = 0; do i = 0, 3; t(i) = sum(der(iy, 0, :)*(4 - i)*dyj**(3 - i)); end do
      der(iy, 1, :) = solve5(M, t)
    end do
    !$omp target enter data map(to: der)
  end subroutine setup_derivatives

  ! Dense 5x5 solve with partial pivoting.
  function solve5(M_in, b_in) result(x)
    real(C_DOUBLE), intent(in) :: M_in(0:4, 0:4), b_in(0:4)
    real(C_DOUBLE) :: x(0:4), M(0:4, 0:4), b(0:4), f, tmp(0:4)
    integer :: i, j, p
    M = M_in; b = b_in
    do i = 0, 4
      p = i - 1 + maxloc(abs(M(i:4, i)), 1)
      if (p /= i) then
        tmp = M(i, :); M(i, :) = M(p, :); M(p, :) = tmp
        f = b(i); b(i) = b(p); b(p) = f
      end if
      do j = i + 1, 4
        f = M(j, i)/M(i, i)
        M(j, i:4) = M(j, i:4) - f*M(i, i:4)
        b(j) = b(j) - f*b(i)
      end do
    end do
    do i = 4, 0, -1
      x(i) = (b(i) - sum(M(i, i + 1:4)*x(i + 1:4)))/M(i, i)
    end do
  end function solve5

  ! The streamwise displacement of the upper image at the current time,
  ! reduced to one box length: exp(-i kx shift) is periodic in it.
  real(C_DOUBLE) function shear_shift_length()
    shear_shift_length = modulo(S*time*ly, lx)
  end function shear_shift_length

  ! Ghost rows of component c of V from its interior rows, at the current time.
  subroutine fill_ghosts(c)
    integer(C_INT), intent(in) :: c
    integer(C_INT) :: ix, iz
    real(C_DOUBLE) :: shift
    complex(C_DOUBLE_COMPLEX) :: ph
    shift = shear_shift_length()
    !$omp target teams distribute parallel do collapse(2) default(none) &
    !$omp shared(V, nx0, nxN, nz, ny, alfa0, shift, c) private(ix, iz, ph)
    do ix = nx0, nxN
      do iz = -nz, nz
        ph = exp(dcmplx(0.0d0, -alfa0*ix*shift))
        V(ny, iz, ix, c) = V(0, iz, ix, c)*ph
        V(ny + 1, iz, ix, c) = V(1, iz, ix, c)*ph
        V(-1, iz, ix, c) = V(ny - 1, iz, ix, c)*conjg(ph)
        V(-2, iz, ix, c) = V(ny - 2, iz, ix, c)*conjg(ph)
      end do
    end do
  end subroutine fill_ghosts

end module hst_derivatives
