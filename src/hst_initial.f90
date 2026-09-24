! The field a run starts from when there is no restart file.
!
! A random, divergence-free field: for each (kx, ky, kz) a random vector
! potential A is drawn and the velocity is its curl, u = i k x A, so the
! divergence vanishes exactly.  y is expanded in Fourier modes ky = 2 pi m/ly
! (at time zero the shear-periodic condition is plain periodicity), and the
! sum over m is evaluated on the y grid.  The amplitude follows
! E(k) = (k/kpeak)^2 exp(1 - (k/kpeak)^2), which peaks at k = kpeak.
!
! The random numbers are keyed on the mode indices (uniform_from_key, from
! channel/src/physics/initial_condition.f90), so the field is a pure
! function of the deck: the same at any rank count.  The ix = 0 plane is
! made Hermitian so that the physical field is real.
module hst_initial

  use, intrinsic :: iso_c_binding
  use hst_params

  implicit none
  private
  public :: generate_initial_field, uniform_from_key

contains

  subroutine generate_initial_field()
    real(C_DOUBLE), parameter :: PI = 3.141592653589793d0
    integer(C_INT) :: ix, iz, iy, m, mmax
    real(C_DOUBLE) :: kx, ky, kz, k, env
    complex(C_DOUBLE_COMPLEX) :: a(3), u(3), ph
    complex(C_DOUBLE_COMPLEX), parameter :: I = (0.0d0, 1.0d0)

    if (has_terminal) write (*, '(A,I0,A,F6.2,A,ES9.2)') &
      '   generating initial field: seed = ', seed, '  kpeak = ', kpeak, '  amplitude = ', amplitude
    V = 0
    mmax = ny/3
    do ix = nx0, nxN
      kx = alfa0*ix
      do iz = -nz, nz
        kz = beta0*iz
        do m = -mmax, mmax
          ky = 2.0d0*PI*m/ly
          k = sqrt(kx*kx + ky*ky + kz*kz)
          if (k == 0.0d0) cycle
          env = (k/kpeak)**2*exp(1.0d0 - (k/kpeak)**2)
          if (env < 1.0d-8) cycle
          call potential(ix, iz, m, a)
          a = amplitude*env*a
          u(1) = I*(ky*a(3) - kz*a(2))
          u(2) = I*(kz*a(1) - kx*a(3))
          u(3) = I*(kx*a(2) - ky*a(1))
          do iy = -2, ny + 1
            ph = exp(I*ky*y(iy))
            V(iy, iz, ix, 1) = V(iy, iz, ix, 1) + u(1)*ph
            V(iy, iz, ix, 2) = V(iy, iz, ix, 2) + u(2)*ph
            V(iy, iz, ix, 3) = V(iy, iz, ix, 3) + u(3)*ph
          end do
        end do
      end do
    end do
  end subroutine generate_initial_field

  ! Random unit-scale complex potential for mode (ix, iz, m).  On the ix = 0
  ! plane the mode (-iz, -m) is the conjugate of (iz, m), which keeps the
  ! physical field real.
  subroutine potential(ix, iz, m, a)
    integer(C_INT), intent(in) :: ix, iz, m
    complex(C_DOUBLE_COMPLEX), intent(out) :: a(3)
    real(C_DOUBLE), parameter :: TWOPI = 6.283185307179586d0
    integer(C_INT) :: c, jz, jm
    real(C_DOUBLE) :: r, phase
    logical :: mirror

    jz = iz; jm = m; mirror = .false.
    if (ix == 0) then
      if (iz < 0 .or. (iz == 0 .and. m < 0)) then
        jz = -iz; jm = -m; mirror = .true.
      end if
    end if
    do c = 1, 3
      r = uniform_from_key(seed, c, jm, jz, ix)
      phase = TWOPI*uniform_from_key(seed, c + 3, jm, jz, ix)
      a(c) = r*exp(dcmplx(0.0d0, phase))
      if (mirror) a(c) = conjg(a(c))
    end do
  end subroutine potential

  ! A uniform deviate in [0, 1) determined by the seed and the mode indices
  ! alone.  An fmix32-shaped avalanche with every intermediate value masked
  ! to 32 bits, so nothing here can overflow (the obvious 64-bit splitmix
  ! finalizer relies on signed overflow and gfortran -O2 miscompiles it).
  pure function uniform_from_key(seed, component, iy, iz, ix) result(u)
    integer(C_INT), intent(in) :: seed, component, iy, iz, ix
    real(C_DOUBLE) :: u
    integer(C_INT64_T) :: h
    integer(C_INT64_T), parameter :: MASK32 = 4294967295_C_INT64_T
    integer(C_INT64_T), parameter :: MIX_A = 2146121005_C_INT64_T
    integer(C_INT64_T), parameter :: MIX_B = 2032289749_C_INT64_T

    h = int(seed, C_INT64_T)
    h = h + 2654435761_C_INT64_T*int(component, C_INT64_T)
    h = h + 2246822519_C_INT64_T*int(iy, C_INT64_T)
    h = h + 3266489917_C_INT64_T*int(iz, C_INT64_T)
    h = h + 668265263_C_INT64_T*int(ix, C_INT64_T)
    h = iand(h, MASK32)
    h = ieor(h, ishft(h, -16))
    h = iand(h*MIX_A, MASK32)
    h = ieor(h, ishft(h, -13))
    h = iand(h*MIX_B, MASK32)
    h = ieor(h, ishft(h, -16))
    u = real(h, C_DOUBLE)*(1.0d0/4294967296.0d0)
  end function uniform_from_key

end module hst_initial
