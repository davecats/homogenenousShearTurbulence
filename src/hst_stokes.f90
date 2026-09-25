! The Stokes layer: an oscillating spanwise mean profile localised at
! mid-box, driven by the body force that makes it an exact solution of the
! mean spanwise momentum equation.  From SLdata.cpl (flags StokesLayer,
! bodyforce, bf_dvw of hst-main).
!
!   W(y, t) = A exp(-s) cos(omega (t - t0) - s),   s = sqrt((eta/delta)^2 + eps^2),
!   eta = y - ly/2,  omega = 2 pi / T,  eps = 0.1  (smoothed |eta|/delta)
!
! is Stokes' second problem about the plane y = ly/2, decaying on both
! sides.  With sl_bodyforce (the CPL default) the mean w equation carries
! the force f = dW/dt - nu d2W/dy2, which for this W is
!
!   f = A exp(-s) [ (2 nu s'^2 - omega) sin(theta) - nu s'' (sin(theta) - cos(theta)) ],
!   theta = omega (t - t0) - s,  s' = (eta/delta^2)/s,  s'' = (eps^2/delta^2)/s^3
!
! (the bodyF of SLdata.cpl in closed form), and the Reynolds-stress
! divergence is dropped from that equation (bf_dvw), so the mean profile is
! exactly W whatever the turbulence does.  The profile is prescribed
! outright during the first two periods after t0 and on the two rows next to
! the box edge (as apply_SL does); without sl_bodyforce it is prescribed at
! every substep instead.  sl_ramp multiplies both by the cubic ramp
! 3 tau^2 - 2 tau^3 over the first period (the CPL smoothStep flag; the
! default is the hard switch-on at t0).
!
! The layer needs the stretched grid (ystretch > 0 in &mesh) to resolve
! delta at mid-box.  Region statistics (inside |eta| < 8 delta and outside)
! go to stokes_runtime.dat, see hst_stats.
!
! Here the force is added D0-weighted, like every other term of the
! equations; SLdata.cpl adds it raw, which is a second-order error in the
! applied force on the D0-weighted equation.
module hst_stokes

  use, intrinsic :: iso_c_binding
  use hst_params

  implicit none
  private
  public :: stokes_active, stokes_setup, stokes_force, stokes_apply, stokes_profile

  real(C_DOUBLE), parameter :: EPS_ABSM = 0.1d0, N_PERIODS = 2.0d0
  real(C_DOUBLE), parameter :: TWOPI = 6.283185307179586d0

contains

  logical function stokes_active()
    stokes_active = (sl_amplitude /= 0.0d0)
  end function stokes_active

  real(C_DOUBLE) function ramp(t)
    real(C_DOUBLE), intent(in) :: t
    real(C_DOUBLE) :: tau
    if (t < sl_start) then
      ramp = 0.0d0
    else if (sl_ramp) then
      tau = min((t - sl_start)/sl_period, 1.0d0)
      ramp = 3.0d0*tau**2 - 2.0d0*tau**3
    else
      ramp = 1.0d0
    end if
  end function ramp

  ! W(y, t), with the ramp.
  real(C_DOUBLE) function stokes_profile(yy, t)
    real(C_DOUBLE), intent(in) :: yy, t
    real(C_DOUBLE) :: s, omega
    omega = TWOPI/sl_period
    s = sqrt(((yy - 0.5d0*ly)/sl_delta)**2 + EPS_ABSM**2)
    stokes_profile = ramp(t)*sl_amplitude*exp(-s)*cos(omega*(t - sl_start) - s)
  end function stokes_profile

  ! f(y, t) = dW/dt - nu d2W/dy2 for the un-ramped W, times the ramp.
  real(C_DOUBLE) function stokes_force_at(yy, t)
    real(C_DOUBLE), intent(in) :: yy, t
    real(C_DOUBLE) :: eta, s, s1, s2, omega, theta
    omega = TWOPI/sl_period
    eta = yy - 0.5d0*ly
    s = sqrt((eta/sl_delta)**2 + EPS_ABSM**2)
    s1 = (eta/sl_delta**2)/s
    s2 = (EPS_ABSM**2/sl_delta**2)/s**3
    theta = omega*(t - sl_start) - s
    stokes_force_at = ramp(t)*sl_amplitude*exp(-s)* &
                      ((2.0d0*ni*s1*s1 - omega)*sin(theta) - ni*s2*(sin(theta) - cos(theta)))
  end function stokes_force_at

  ! Marks the rows inside the layer (|eta| < 8 delta) for the statistics.
  subroutine stokes_setup()
    integer(C_INT) :: iy
    inlayer = 0
    if (.not. stokes_active()) return
    do iy = 0, ny - 1
      if (abs(y(iy) - 0.5d0*ly) < 8.0d0*sl_delta) inlayer(iy) = 1
    end do
    !$omp target update to(inlayer)
    if (has_terminal) write (*, '(A,I0,A,I0,A)') '   Stokes layer: ', sum(inlayer), ' of ', ny, ' rows inside |y - ly/2| < 8 delta'
  end subroutine stokes_setup

  ! The body force on the mean w equation at time t, with periodic ghost
  ! rows (it is applied D0-weighted), on the device.
  subroutine stokes_force(t)
    real(C_DOUBLE), intent(in) :: t
    integer(C_INT) :: iy
    if (.not. (stokes_active() .and. sl_bodyforce)) return
    do iy = 0, ny - 1
      fy(iy) = stokes_force_at(y(iy), t)
    end do
    fy(-2) = fy(ny - 2); fy(-1) = fy(ny - 1); fy(ny) = fy(0); fy(ny + 1) = fy(1)
    !$omp target update to(fy)
  end subroutine stokes_force

  ! Prescribes the mean spanwise profile (the (0,0) mode of w, real) on the
  ! rank that owns it: everywhere without the body force or during the first
  ! two periods, and on the two rows next to the box edge always.  Leaves
  ! the device copy of that mode updated; the caller refills the ghosts.
  subroutine stokes_apply()
    integer(C_INT) :: iy
    logical :: whole
    if (.not. stokes_active() .or. .not. has_average) return
    whole = (.not. sl_bodyforce) .or. (time - sl_start <= N_PERIODS*sl_period)
    !$omp target update from(V(:, 0, 0, 3))
    do iy = 0, ny - 1
      if (whole .or. iy == 0 .or. iy == ny - 1) V(iy, 0, 0, 3) = dcmplx(stokes_profile(y(iy), time), 0.0d0)
    end do
    !$omp target update to(V(:, 0, 0, 3))
  end subroutine stokes_apply

end module hst_stokes
