! The equations.
!
! Fluctuations about U = S*y in velocity-vorticity form.  Per Fourier mode
! (alfa, beta), with D the compact operators of hst_derivatives and
! k2 = alfa^2 + beta^2, the unknowns are
!   d2v = (D2 - k2 D0) v                       (D0-weighted Laplacian of v)
!   eta = D0 (i beta u - i alfa w)             (D0-weighted vertical vorticity)
! and for the (0,0) mode eta packs the two real mean profiles (u, w).
!
!   d/dt d2v = ni (D4 - 2 k2 D2 + k2^2 D0) v  + nonlinear terms
!   d/dt eta = ni (D2 - k2 D0) eta            + nonlinear terms - S i beta D0 v
!
! The mean-shear advection S y d/dx is not in these: it is integrated
! exactly by the phase exp(-i alfa S y dt) applied after each substep
! (shear_shift).  Time stepping: RK3 (Rai-Moin) for the explicit terms,
! Crank-Nicolson for the viscous ones, coefficients ODE = RK_rai(:, i):
!   rhs = ODE(1) unkn/deltat + impl + ODE(2) expl_new - ODE(3) expl_old
! then (ODE(1)/deltat - viscous) unkn_new = rhs.
!
! From channel/src/physics/channel_equations.fypp with the fypp macros
! written out, the flow-rate correction and the scalars removed, and the
! two shear terms added.  The nonlinear terms follow S1data.cpl.
module hst_equations

  use, intrinsic :: iso_c_binding
  use hst_params
  use hst_derivatives, only: fill_ghosts
  use hst_linsolve, only: solve_component, apply_dy, KIND_D2V, KIND_ETA
  use hst_fft, only: VVdz
  use hst_transforms, only: transform_to_physical, build_products, products_to_spectral, compute_cfl

  implicit none
  private
  public :: timestep, buildrhs_prepare, buildrhs, shear_shift, linsolve

contains

  ! One full time step of length deltat: three RK substeps, each
  !   transform V to physical space (ghost rows included)
  !   build the right-hand sides from V and from the six products
  !   shift them to the new time frame (exact mean-shear advection)
  !   advance time; solve the implicit systems; recover u, w; fill ghosts
  ! The CFL number of the third substep is left in cfl (rank-local).
  ! This is the loop of scddns.cpl.
  subroutine timestep()
    integer(C_INT) :: i, m
    real(C_DOUBLE) :: dt_sub
    do i = 1, 3
      if (.not. linear .or. i == 3) call transform_to_physical()
      if (i == 3) call compute_cfl()
      call buildrhs_prepare(RK_rai(:, i))
      if (.not. linear) then
        do m = 1, 6
          call build_products(m)
          call products_to_spectral()
          call buildrhs(RK_rai(:, i), m)
        end do
      end if
      dt_sub = 2.0d0/RK_rai(1, i)*deltat
      call shear_shift(dt_sub)
      time = time + dt_sub
      call linsolve(RK_rai(1, i)/deltat)
    end do
  end subroutine timestep

  ! The parts of the right-hand side that depend on V only: the unknowns
  ! divided by deltat, the explicit half of the viscous terms, the carried
  ! explicit term of the previous substep, and the mean-shear tilting term.
  ! Writes memrhs, then copies it into V(:, :, :, 1:2), which from here on
  ! hold the right-hand sides until linsolve overwrites them.
  subroutine buildrhs_prepare(ODE)
    real(C_DOUBLE), intent(in) :: ODE(3)
    integer(C_INT) :: ix, iy, iz, j
    complex(C_DOUBLE_COMPLEX) :: unkn, impl, expl, f
    real(C_DOUBLE) :: helm, biharm

    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(V, memrhs, oldrhs, der, k2, ialfa, ibeta, ni, S, deltat, ODE, nx0, nxN, nz, ny) &
    !$omp private(ix, iy, iz, j, unkn, impl, expl, f, helm, biharm)
    do ix = nx0, nxN
      do iz = -nz, nz
        do iy = 0, ny - 1
          ! v equation
          unkn = 0.0d0; impl = 0.0d0
          do j = -2, 2
            helm = der(iy, 2, j) - k2(iz, ix)*der(iy, 0, j)
            biharm = der(iy, 3, j) - 2.0d0*k2(iz, ix)*der(iy, 2, j) + k2(iz, ix)*k2(iz, ix)*der(iy, 0, j)
            unkn = unkn + helm*V(iy + j, iz, ix, 2)
            impl = impl + ni*biharm*V(iy + j, iz, ix, 2)
          end do
          memrhs(iy, iz, ix, 2) = ODE(1)*unkn/deltat + impl - ODE(3)*oldrhs(iy, iz, ix, 2)
          oldrhs(iy, iz, ix, 2) = 0.0d0
          ! eta equation (mean mode: u and w packed as one complex number)
          unkn = 0.0d0; impl = 0.0d0; expl = 0.0d0
          do j = -2, 2
            helm = der(iy, 2, j) - k2(iz, ix)*der(iy, 0, j)
            if (ix == 0 .and. iz == 0) then
              f = dcmplx(dreal(V(iy + j, iz, ix, 1)), dreal(V(iy + j, iz, ix, 3)))
              unkn = unkn + der(iy, 0, j)*f
              impl = impl + ni*der(iy, 2, j)*f
            else
              f = ibeta(iz)*V(iy + j, iz, ix, 1) - ialfa(ix)*V(iy + j, iz, ix, 3)
              unkn = unkn + der(iy, 0, j)*f
              impl = impl + ni*helm*f
            end if
            expl = expl - S*ibeta(iz)*der(iy, 0, j)*V(iy + j, iz, ix, 2)   ! tilting of the mean vorticity
          end do
          memrhs(iy, iz, ix, 1) = ODE(1)*unkn/deltat + impl - ODE(3)*oldrhs(iy, iz, ix, 1) + ODE(2)*expl
          oldrhs(iy, iz, ix, 1) = expl
        end do
      end do
    end do
    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(V, memrhs, nx0, nxN, nz, ny) private(ix, iy, iz)
    do ix = nx0, nxN
      do iz = -nz, nz
        do iy = 0, ny - 1
          V(iy, iz, ix, 1) = memrhs(iy, iz, ix, 1)
          V(iy, iz, ix, 2) = memrhs(iy, iz, ix, 2)
        end do
      end do
    end do
  end subroutine buildrhs_prepare

  ! Nonlinear contribution of product m (1..6 = uu, vv, ww, uv, vw, uw),
  ! which products_to_spectral has left in VVdz, to the two right-hand
  ! sides.  With -div(u u) written per component,
  !   rhsu = -i alfa uu - D1 uv - i beta uw
  !   rhsw = -i alfa uw - D1 vw - i beta ww
  !   rhsv = -i alfa uv - D1 vv - i beta vw
  ! the eta equation gets i beta rhsu - i alfa rhsw and the d2v equation
  ! D1 (i alfa rhsu + i beta rhsw) - k2 rhsv, split here by product.
  subroutine buildrhs(ODE, m)
    real(C_DOUBLE), intent(in) :: ODE(3)
    integer(C_INT), intent(in) :: m
    integer(C_INT) :: ix, iy, iz, j
    complex(C_DOUBLE_COMPLEX) :: d0, d1, d2, rhsu, rhsw, expl, e

    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(V, oldrhs, VVdz, der, izd, k2, ialfa, ibeta, ODE, m, nx0, nxN, nz, ny) &
    !$omp private(ix, iy, iz, j, d0, d1, d2, rhsu, rhsw, expl, e)
    do ix = nx0, nxN
      do iz = -nz, nz
        do iy = 0, ny - 1
          d0 = 0.0d0; d1 = 0.0d0; d2 = 0.0d0
          do j = -2, 2
            d0 = d0 + der(iy, 0, j)*VVdz(izd(iz) + 1, ix - nx0 + 1, iy + j)
            d1 = d1 + der(iy, 1, j)*VVdz(izd(iz) + 1, ix - nx0 + 1, iy + j)
            d2 = d2 + der(iy, 2, j)*VVdz(izd(iz) + 1, ix - nx0 + 1, iy + j)
          end do
          select case (m)
          case (1)   ! uu
            rhsu = -ialfa(ix)*d0; rhsw = 0.0d0
            expl = ialfa(ix)*ialfa(ix)*d1
          case (2)   ! vv
            rhsu = 0.0d0; rhsw = 0.0d0
            expl = k2(iz, ix)*d1
          case (3)   ! ww
            rhsu = 0.0d0; rhsw = -ibeta(iz)*d0
            expl = ibeta(iz)*ibeta(iz)*d1
          case (4)   ! uv
            rhsu = -d1; rhsw = 0.0d0
            expl = ialfa(ix)*d2 + ialfa(ix)*k2(iz, ix)*d0
          case (5)   ! vw
            rhsu = 0.0d0; rhsw = -d1
            expl = ibeta(iz)*d2 + ibeta(iz)*k2(iz, ix)*d0
          case default   ! uw
            rhsu = -ibeta(iz)*d0; rhsw = -ialfa(ix)*d0
            expl = 2.0d0*ialfa(ix)*ibeta(iz)*d1
          end select
          V(iy, iz, ix, 2) = V(iy, iz, ix, 2) + ODE(2)*expl
          oldrhs(iy, iz, ix, 2) = oldrhs(iy, iz, ix, 2) + expl
          if (ix == 0 .and. iz == 0) then
            e = dcmplx(dreal(rhsu), dreal(rhsw))
          else
            e = ibeta(iz)*rhsu - ialfa(ix)*rhsw
          end if
          V(iy, iz, ix, 1) = V(iy, iz, ix, 1) + ODE(2)*e
          oldrhs(iy, iz, ix, 1) = oldrhs(iy, iz, ix, 1) + e
        end do
      end do
    end do
  end subroutine buildrhs

  ! Exact integration of the mean-shear advection over a substep of length
  ! dt_sub: the right-hand sides (in V(:, :, :, 1:2)) and the carried
  ! explicit terms move to the new time frame.
  subroutine shear_shift(dt_sub)
    real(C_DOUBLE), intent(in) :: dt_sub
    integer(C_INT) :: ix, iy, iz
    complex(C_DOUBLE_COMPLEX) :: f
    ! (no default(none): nvfortran 25.9 rejects the grid array y in a
    !  shared clause here, although it accepts it elsewhere)
    !$omp target teams distribute parallel do collapse(3) &
    !$omp shared(V, oldrhs, y, alfa0, S, dt_sub, nx0, nxN, nz, ny) private(ix, iy, iz, f)
    do ix = nx0, nxN
      do iz = -nz, nz
        do iy = 0, ny - 1
          f = exp(dcmplx(0.0d0, -alfa0*ix*S*y(iy)*dt_sub))
          V(iy, iz, ix, 1) = V(iy, iz, ix, 1)*f
          V(iy, iz, ix, 2) = V(iy, iz, ix, 2)*f
          oldrhs(iy, iz, ix, 1) = oldrhs(iy, iz, ix, 1)*f
          oldrhs(iy, iz, ix, 2) = oldrhs(iy, iz, ix, 2)*f
        end do
      end do
    end do
  end subroutine shear_shift

  ! Implicit solves, then u and w from continuity and the definition of eta:
  !   i alfa u + dv/dy + i beta w = 0,   i beta u - i alfa w = eta
  ! Leaves V with all ghost rows filled at the current time.
  subroutine linsolve(lambda)
    real(C_DOUBLE), intent(in) :: lambda
    integer(C_INT) :: ix, iy, iz
    complex(C_DOUBLE_COMPLEX) :: temp

    call solve_component(KIND_D2V, lambda, V(:, :, :, 2))
    call solve_component(KIND_ETA, lambda, V(:, :, :, 1))
    call fill_ghosts(2)
    call apply_dy(2, V(:, :, :, 3))            ! V(:, :, :, 3) = dv/dy
    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(V, k2, ialfa, ibeta, nx0, nxN, nz, ny) private(ix, iy, iz, temp)
    do ix = nx0, nxN
      do iz = -nz, nz
        do iy = 0, ny - 1
          if (ix == 0 .and. iz == 0) then
            temp = V(iy, iz, ix, 1)
            V(iy, iz, ix, 1) = dcmplx(dreal(temp), 0.0d0)
            V(iy, iz, ix, 3) = dcmplx(dimag(temp), 0.0d0)
          else
            temp = (ialfa(ix)*V(iy, iz, ix, 3) - ibeta(iz)*V(iy, iz, ix, 1))/k2(iz, ix)
            V(iy, iz, ix, 3) = (ibeta(iz)*V(iy, iz, ix, 3) + ialfa(ix)*V(iy, iz, ix, 1))/k2(iz, ix)
            V(iy, iz, ix, 1) = temp
          end if
        end do
      end do
    end do
    call fill_ghosts(1)
    call fill_ghosts(3)
  end subroutine linsolve

end module hst_equations
