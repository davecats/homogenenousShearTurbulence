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
  public :: setup_derivatives, fill_ghosts, fill_ghosts_field, shear_shifts
  public :: s2_of, s2_integral, gamma_y_of

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

  !------------------------------------------------- unsteady spanwise shear ----
  ! S2(t) as in S2data.cpl: zero before s2_start, then s2_amplitude times
  ! sin(omega (t - s2_start)) with omega = 2 pi / s2_period, or constant
  ! when s2_period = 0.
  real(C_DOUBLE) function s2_of(t)
    real(C_DOUBLE), intent(in) :: t
    real(C_DOUBLE), parameter :: TWOPI = 6.283185307179586d0
    s2_of = 0.0d0
    if (s2_amplitude == 0.0d0 .or. t < s2_start) return
    if (s2_period > 0.0d0) then
      s2_of = s2_amplitude*sin(TWOPI/s2_period*(t - s2_start))
    else
      s2_of = s2_amplitude
    end if
  end function s2_of

  ! Integral of S2 from t1 to t2 (t1 <= t2), closed form.
  real(C_DOUBLE) function s2_integral(t1, t2)
    real(C_DOUBLE), intent(in) :: t1, t2
    real(C_DOUBLE), parameter :: TWOPI = 6.283185307179586d0
    real(C_DOUBLE) :: a, b, omega
    s2_integral = 0.0d0
    if (s2_amplitude == 0.0d0 .or. t2 <= s2_start) return
    a = max(t1, s2_start); b = t2
    if (s2_period > 0.0d0) then
      omega = TWOPI/s2_period
      s2_integral = s2_amplitude/omega*(cos(omega*(a - s2_start)) - cos(omega*(b - s2_start)))
    else
      s2_integral = s2_amplitude*(b - a)
    end if
  end function s2_integral

  ! The spanwise displacement of the upper image, gamma_y(t) = int S2 dt.
  real(C_DOUBLE) function gamma_y_of(t)
    real(C_DOUBLE), intent(in) :: t
    gamma_y_of = s2_integral(s2_start, t)
  end function gamma_y_of

  ! The streamwise and spanwise displacements of the upper image at time t,
  ! each reduced to one box length (the wrap phase is periodic in them):
  !   f_hat(ny) = f_hat(0) * exp(-i (kx shift_x + kz shift_z))
  subroutine shear_shifts(t, shift_x, shift_z)
    real(C_DOUBLE), intent(in) :: t
    real(C_DOUBLE), intent(out) :: shift_x, shift_z
    shift_x = modulo(S*t*ly, lx)
    shift_z = modulo(gamma_y_of(t)*ly, lz)
  end subroutine shear_shifts

  ! Ghost rows of component c of V from its interior rows, at the current time.
  subroutine fill_ghosts(c)
    integer(C_INT), intent(in) :: c
    real(C_DOUBLE) :: shift_x, shift_z
    call shear_shifts(time, shift_x, shift_z)
    call fill_ghosts_field(V(:, :, :, c), shift_x, shift_z)
  end subroutine fill_ghosts

  ! Ghost rows of any field with the layout of a component of V, for the
  ! displacements shift_x, shift_z of the upper image.
  subroutine fill_ghosts_field(field, shift_x, shift_z)
    complex(C_DOUBLE_COMPLEX), intent(inout) :: field(ny0 - 2:, -nz:, nx0:)
    real(C_DOUBLE), intent(in) :: shift_x, shift_z
    integer(C_INT) :: ix, iz
    complex(C_DOUBLE_COMPLEX) :: ph
    !$omp target teams distribute parallel do collapse(2) default(none) &
    !$omp shared(field, nx0, nxN, nz, ny, alfa0, beta0, shift_x, shift_z) private(ix, iz, ph)
    do ix = nx0, nxN
      do iz = -nz, nz
        ph = exp(dcmplx(0.0d0, -(alfa0*ix*shift_x + beta0*iz*shift_z)))
        field(ny, iz, ix) = field(0, iz, ix)*ph
        field(ny + 1, iz, ix) = field(1, iz, ix)*ph
        field(-1, iz, ix) = field(ny - 1, iz, ix)*conjg(ph)
        field(-2, iz, ix) = field(ny - 2, iz, ix)*conjg(ph)
      end do
    end do
  end subroutine fill_ghosts_field

end module hst_derivatives
