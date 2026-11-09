! Large-scale damping (DAMPING.md).
!
! A body force -g(kx, kz) u on the Fourier columns inside a region of the
! x-z wavenumber plane defined in the shear frame: with e_par the direction
! of the shear vector (S, S2(t)) (or a fixed angle), k_par = k . e_par,
! k_perp = k . e_perp and
!
!   r = ( (|k_par|/damp_kpar)^p + (|k_perp|/damp_kperp)^p )^(1/p),   p = damp_norm,
!
! the columns with r < 1 are damped (r = 2: an ellipse with the semi-axes
! damp_kpar, damp_kperp; p large: the rectangle).  Modes on the boundary
! (r = 1 within 1e-6) are kept, so with the semi-axes 2 pi/Lx, 2 pi/Lz of a
! small box the small box's own modes are the first ones kept.
!
!   soft (damp_rate > 0): g = damp_rate inside, falling to zero through a
!     smoothstep over the last damp_edge of r (damp_edge = 0: a step);
!     treated explicitly with the other explicit terms (RK3), so
!     damp_rate * deltat must stay below about 1;
!   hard (damp_hard): the columns are zeroed every substep (their right-hand
!     sides are dropped before the implicit solve), the limit damp_rate ->
!     infinity: a Galerkin truncation, which removes no power.
!
! The (0,0) column (the mean profiles) is damped only with damp_mean.  With
! damp_follow the region turns with the shear direction atan2(S2(t), S)
! and gdamp is rebuilt every substep in which the direction changed;
! otherwise it is built once.
!
! Diagnostics (hst_stats, damping_runtime.dat): the energy in the damped
! columns and the power 2 g <u.u> they lose.
module hst_damping

  use, intrinsic :: iso_c_binding
  use hst_params
  use hst_derivatives, only: s2_of

  implicit none
  private
  public :: damping_active, damping_setup, damping_update, damping_theta, damping_mask, damping_soft_factor, ncol_damped

  real(C_DOUBLE), parameter :: PI = 3.141592653589793d0
  real(C_DOUBLE), save :: theta_last = -huge(1.0d0)
  integer(C_INT), save :: ncol_damped = 0        ! this rank's damped columns

contains

  logical function damping_active()
    damping_active = (damp_kpar > 0.0d0 .and. damp_kperp > 0.0d0 .and. (damp_hard .or. damp_rate > 0.0d0))
  end function damping_active

  ! 1 when the explicit damping term is in use, 0 for the hard mask or no damping.
  real(C_DOUBLE) function damping_soft_factor()
    damping_soft_factor = 0.0d0
    if (damping_active() .and. .not. damp_hard) damping_soft_factor = 1.0d0
  end function damping_soft_factor

  ! Direction of e_par (radians from x towards z) at time t.
  real(C_DOUBLE) function damping_theta(t)
    real(C_DOUBLE), intent(in) :: t
    if (damp_follow) then
      damping_theta = atan2(s2_of(t), S)
    else
      damping_theta = damp_angle*PI/180.0d0
    end if
  end function damping_theta

  ! Allocates gdamp (always: the equations read it) and builds it for the
  ! start time.
  subroutine damping_setup()
    allocate (gdamp(-nz:nz, nx0:nxN))
    gdamp = 0.0d0
    !$omp target enter data map(to: gdamp)
    if (damping_active()) call damping_update(time)
    if (has_terminal .and. damping_active()) write (*, '(A,F8.3,A,L1,A,F8.3,A,I0)') &
      '   damping: direction at start ', damping_theta(time)*180.0d0/PI, ' deg,  hard = ', damp_hard, &
      '  rate = ', damp_rate, '  damped columns on rank 0: ', ncol_damped
  end subroutine damping_setup

  ! gdamp for the shear direction at time t; nothing to do when the
  ! direction has not changed since the last call.
  subroutine damping_update(t)
    real(C_DOUBLE), intent(in) :: t
    integer(C_INT) :: ix, iz, n
    real(C_DOUBLE) :: theta, c, s, kx, kz, kpar, kperp, r, g, x

    theta = damping_theta(t)
    if (theta == theta_last) return
    theta_last = theta
    c = cos(theta); s = sin(theta)
    n = 0
    do ix = nx0, nxN
      kx = alfa0*ix
      do iz = -nz, nz
        kz = beta0*iz
        kpar = kx*c + kz*s
        kperp = -kx*s + kz*c
        r = ((abs(kpar)/damp_kpar)**damp_norm + (abs(kperp)/damp_kperp)**damp_norm)**(1.0d0/damp_norm)
        g = 0.0d0
        if (r < 1.0d0 - 1.0d-6) then
          if (damp_hard) then
            g = 1.0d0
          else if (r <= 1.0d0 - damp_edge) then
            g = damp_rate
          else
            x = (1.0d0 - r)/damp_edge
            g = damp_rate*x*x*(3.0d0 - 2.0d0*x)
          end if
        end if
        if (ix == 0 .and. iz == 0 .and. .not. damp_mean) g = 0.0d0
        gdamp(iz, ix) = g
        if (g > 0.0d0) n = n + 1
      end do
    end do
    ncol_damped = n
    !$omp target update to(gdamp)
  end subroutine damping_update

  ! Hard mask: drops the right-hand sides of the masked columns before the
  ! implicit solve, so the solve returns zero for them.
  subroutine damping_mask()
    integer(C_INT) :: ix, iy, iz
    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(rhs, oldrhs, gdamp, nx0, nxN, nz, ny0, nyN) private(ix, iy, iz)
    do ix = nx0, nxN
      do iz = -nz, nz
        do iy = ny0, nyN
          if (gdamp(iz, ix) > 0.0d0) then
            rhs(iy, iz, ix, 1) = 0.0d0; rhs(iy, iz, ix, 2) = 0.0d0
            oldrhs(iy, iz, ix, 1) = 0.0d0; oldrhs(iy, iz, ix, 2) = 0.0d0
          end if
        end do
      end do
    end do
  end subroutine damping_mask

end module hst_damping
