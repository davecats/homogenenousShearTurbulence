! The spectral <-> physical round trip.
!
! transform_to_physical: V -> rVVdx(:, :, :, 1:3), i.e. u, v, w on the
!   dealiased physical grid, all rows -2..ny+1 (the ghost rows included, so
!   that products in the ghost rows are the shear-periodic images of the
!   products, which is what the stencils near the box edges need).
! build_products(m) + products_to_spectral: one product of two components
!   back to spectral space, left in VVdz for the caller to accumulate.
! compute_cfl: read off rVVdx while it exists.
!
! From channel/src/numerics/channel_transforms.f90 without the overlapped
! double buffering.
module hst_transforms

  use, intrinsic :: iso_c_binding
  use hst_params
  use hst_mpi, only: transpose_zTOx, transpose_xTOz
  use hst_fft

  implicit none
  private
  public :: transform_to_physical, build_products, products_to_spectral, vvdz_to_field, compute_cfl

contains

  ! V(:, :, :, m) into the z-padded pencil: modes 0..nz first, then zeros,
  ! then -nz..-1 at the end (izd does this index map).
  subroutine assemble_vvdz(m)
    integer(C_INT), intent(in) :: m
    integer(C_INT) :: i, j, k, y_first, y_last
    y_first = ny0 - 2
    y_last = nyN + 2
    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(V, VVdz, nx0, nxN, nz, nzd, y_first, y_last, m) private(i, j, k)
    do i = y_first, y_last
      do j = nx0, nxN
        do k = 1, nzd
          if (k <= nz + 1) then
            VVdz(k, j - nx0 + 1, i) = V(i, k - 1, j, m)
          else if (k <= nzd - nz) then
            VVdz(k, j - nx0 + 1, i) = 0.0d0
          else
            VVdz(k, j - nx0 + 1, i) = V(i, k - nzd - 1, j, m)
          end if
        end do
      end do
    end do
  end subroutine assemble_vvdz

  ! The x modes above nx are padding for the 3/2 rule: zero them before the
  ! real transform.
  subroutine zero_vvdx_padding()
    integer(C_INT) :: i, j, k, y_first, y_last
    y_first = ny0 - 2
    y_last = nyN + 2
    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(VVdx, nx, nxd, nzB, y_first, y_last) private(i, j, k)
    do i = y_first, y_last
      do j = 1, nzB
        do k = nx + 2, nxd + 1
          VVdx(k, j, i) = 0.0d0
        end do
      end do
    end do
  end subroutine zero_vvdx_padding

  subroutine transform_to_physical()
    integer(C_INT) :: m
    do m = 1, 3
      call assemble_vvdz(m)
      call IFT(VVdz)
      call transpose_zTOx(VVdz, VVdx)
      call zero_vvdx_padding()
      call RFT(VVdx, rVVdx(:, :, :, m))
    end do
  end subroutine transform_to_physical

  ! products = a*b*factor with (a, b) = uu, vv, ww, uv, vw, uw for m = 1..6.
  ! factor = 1/(2 nxd nzd) is the normalisation of the inverse transforms.
  subroutine build_products(m)
    integer(C_INT), intent(in) :: m
    integer(C_INT) :: i, j, k, a, b, y_first, y_last
    integer(C_INT), parameter :: first(6) = [1, 2, 3, 1, 2, 1], second(6) = [1, 2, 3, 2, 3, 3]
    y_first = ny0 - 2
    y_last = nyN + 2
    a = first(m)
    b = second(m)
    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(rVVdx, products, nxd, nzB, y_first, y_last, factor, a, b) private(i, j, k)
    do i = y_first, y_last
      do j = 1, nzB
        do k = 1, 2*nxd
          products(k, j, i) = rVVdx(k, j, i, a)*rVVdx(k, j, i, b)*factor
        end do
      end do
    end do
  end subroutine build_products

  ! products (physical) -> VVdz (spectral z-pencil).  The result for mode
  ! (iz, ix) is VVdz(izd(iz) + 1, ix - nx0 + 1, iy).
  subroutine products_to_spectral()
    call HFT(products, VVdx)
    call transpose_xTOz(VVdx, VVdz)
    call FFT(VVdz)
  end subroutine products_to_spectral

  ! Unpack VVdz into a field with the layout of V(:, :, :, 1).
  subroutine vvdz_to_field(field)
    complex(C_DOUBLE_COMPLEX), intent(out) :: field(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
    integer(C_INT) :: ix, iy, iz, y_first, y_last
    y_first = ny0 - 2
    y_last = nyN + 2
    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(VVdz, izd, field, nx0, nxN, nz, y_first, y_last) private(ix, iy, iz)
    do ix = nx0, nxN
      do iy = y_first, y_last
        do iz = -nz, nz
          field(iy, iz, ix) = VVdz(izd(iz) + 1, ix - nx0 + 1, iy)
        end do
      end do
    end do
  end subroutine vvdz_to_field

  ! Rank-local CFL estimate from the physical velocities; the caller reduces
  ! it over ranks.  Only the fluctuations count: the mean shear advection is
  ! integrated analytically and imposes no step restriction.
  subroutine compute_cfl()
    integer(C_INT) :: i, j, k, y_first, y_last
    real(C_DOUBLE) :: tmp
    y_first = ny0
    y_last = nyN
    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(rVVdx, dx, dyl, dz, y_first, y_last, nxd, nzB) private(i, j, k, tmp) reduction(max:cfl)
    do i = y_first, y_last
      do k = 1, nzB
        do j = 1, 2*nxd
          tmp = abs(rVVdx(j, k, i, 1))/dx + abs(rVVdx(j, k, i, 2))/dyl(i) + abs(rVVdx(j, k, i, 3))/dz
          cfl = max(cfl, tmp)
        end do
      end do
    end do
  end subroutine compute_cfl

end module hst_transforms
