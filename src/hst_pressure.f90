! Pressure, computed online at snapshot times.
!
! The fluctuation pressure in the mean shear U = S y satisfies, per mode,
!   (D2 - k2 D0) p = D0 [ alfa^2 uu + beta^2 ww + 2 alfa beta uw - 2 (S i alfa + S2 i beta) v ]
!                    - D2 vv - 2 i alfa D1 uv - 2 i beta D1 vw
! (the D0-weighted form of  lap p = -d_i d_j (u_i u_j) - 2 S dv/dx), with
! the same shear-periodic wrap as the velocity, solved by the cyclic line
! solver.  The (0,0) mode is singular: p_00 = -<vv>_xz + const, with the
! constant fixed by zero mean.
!
! The six products are the ones the solver forms every substep; here they
! are formed once more on the freshly solved field.  The result goes to
! Dati.cart.<i>.p.out next to the velocity snapshot.  Reduced from
! channel/src/io/pressure_output.fypp: no wall rows, no dp/dy file.
module hst_pressure

  use, intrinsic :: iso_c_binding
  use mpi_f08
  use hst_params
  use hst_fft, only: VVdz
  use hst_transforms, only: transform_to_physical, build_products, products_to_spectral
  use hst_linsolve, only: solve_component, KIND_POISSON
  use hst_io, only: field_write
  use hst_derivatives, only: s2_of

  implicit none
  private
  public :: compute_pressure, write_pressure

contains

  ! p (layout of one component of V, on the device) from the current V,
  ! whose ghost rows must be filled.  Uses the products buffers.
  subroutine compute_pressure(p)
    complex(C_DOUBLE_COMPLEX), intent(inout) :: p(ny0 - 2:, -nz:, nx0:)
    integer(C_INT) :: m, ix, iy, iz, j
    complex(C_DOUBLE_COMPLEX) :: d0, d1, d2, term
    real(C_DOUBLE) :: pmean, pmean_global, s2now
    integer :: ierr

    s2now = s2_of(time)
    call transform_to_physical()
    do m = 1, 6
      call build_products(m)
      call products_to_spectral()
      !$omp target teams distribute parallel do collapse(3) default(none) &
      !$omp shared(p, VVdz, V, der, izd, ialfa, ibeta, S, s2now, m, nx0, nxN, nz, ny) &
      !$omp private(ix, iy, iz, j, d0, d1, d2, term)
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
              term = -ialfa(ix)*ialfa(ix)*d0
            case (2)   ! vv
              term = -d2
            case (3)   ! ww
              term = -ibeta(iz)*ibeta(iz)*d0
            case (4)   ! uv
              term = -2.0d0*ialfa(ix)*d1
            case (5)   ! vw
              term = -2.0d0*ibeta(iz)*d1
            case default   ! uw
              term = -2.0d0*ialfa(ix)*ibeta(iz)*d0
            end select
            if (m == 1) then
              ! first product: start the sum with the mean-shear term
              d0 = 0.0d0
              do j = -2, 2
                d0 = d0 + der(iy, 0, j)*V(iy + j, iz, ix, 2)
              end do
              p(iy, iz, ix) = term - 2.0d0*(S*ialfa(ix) + s2now*ibeta(iz))*d0
            else
              p(iy, iz, ix) = p(iy, iz, ix) + term
            end if
            ! mean mode: p_00 = -<vv>, the raw product, not its stencil
            if (ix == 0 .and. iz == 0 .and. m == 2) p(iy, iz, ix) = -dreal(VVdz(1, 1, iy))
          end do
        end do
      end do
    end do
    call solve_component(KIND_POISSON, 0.0d0, p)
    ! zero-mean gauge for the (0,0) mode
    pmean = 0.0d0
    if (has_average) then
      !$omp target teams distribute parallel do default(none) shared(p, ny) private(iy) reduction(+:pmean)
      do iy = 0, ny - 1
        pmean = pmean + dreal(p(iy, 0, 0))
      end do
      pmean = pmean/ny
      !$omp target teams distribute parallel do default(none) shared(p, ny, pmean) private(iy)
      do iy = 0, ny - 1
        p(iy, 0, 0) = dcmplx(dreal(p(iy, 0, 0)) - pmean, 0.0d0)
      end do
    end if
  end subroutine compute_pressure

  ! Pressure of the current field to a file.  Uses memrhs(:, :, :, 2) as
  ! scratch (free outside a substep).
  subroutine write_pressure(filename)
    character(len=*), intent(in) :: filename
    if (has_terminal) print '(A,F12.5)', '   writing '//trim(filename)//' at time', time
    call compute_pressure(memrhs(:, :, :, 2))
    !$omp target update from(memrhs)
    call field_write(filename, memrhs(:, :, :, 2))
  end subroutine write_pressure

end module hst_pressure
