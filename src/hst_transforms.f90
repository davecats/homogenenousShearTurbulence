! The spectral <-> physical round trip.
!
! transform_to_physical: V -> rVVdx(:, :, :, 1:3), i.e. u, v, w on the
!   dealiased physical grid, all rows -2..ny+1 (the ghost rows included, so
!   that products in the ghost rows are the shear-periodic images of the
!   products, which is what the stencils near the box edges need).
! build_products(g) + products_to_spectral: three of the six products
!   (group 1: uu, vv, ww; group 2: uv, vw, uw) back to spectral space, left
!   in VVdz(:, :, :, 1:3) for the caller to accumulate.
! compute_cfl: read off rVVdx while it exists.
!
! The transforms and transposes go field by field, so that the alltoall of
! one field (started by transpose_*_start, waited for by *_finish, two in
! flight at most: hst_mpi) overlaps the transforms of its neighbours:
!   z transform 1, start 1, z transform 2, start 2, finish 1, x transform 1,
!   z transform 3, start 3, finish 2, x transform 2, finish 3, x transform 3.
! This is the double buffering of channel/src/numerics/channel_transforms.f90.
module hst_transforms

  use, intrinsic :: iso_c_binding
  use hst_params
  use hst_mpi, only: transpose_zTOx_start, transpose_zTOx_finish, transpose_xTOz_start, transpose_xTOz_finish
  use hst_fft

  implicit none
  private
  public :: transform_to_physical, build_products, products_to_spectral, vvdz_to_field, compute_cfl

contains

  ! V into the z-padded pencil: modes 0..nz first, then zeros, then -nz..-1
  ! at the end (izd does this index map).
  subroutine assemble_vvdz()
    integer(C_INT) :: i, j, k, m, y_first, y_last
    y_first = ny0 - 2
    y_last = nyN + 2
    !$omp target teams distribute parallel do collapse(4) default(none) &
    !$omp shared(V, VVdz, nx0, nxN, nz, nzd, y_first, y_last) private(i, j, k, m)
    do m = 1, 3
      do i = y_first, y_last
        do j = nx0, nxN
          do k = 1, nzd
            if (k <= nz + 1) then
              VVdz(k, j - nx0 + 1, i, m) = V(i, k - 1, j, m)
            else if (k <= nzd - nz) then
              VVdz(k, j - nx0 + 1, i, m) = 0.0d0
            else
              VVdz(k, j - nx0 + 1, i, m) = V(i, k - nzd - 1, j, m)
            end if
          end do
        end do
      end do
    end do
  end subroutine assemble_vvdz

  ! The x modes above nx of field m are padding for the 3/2 rule: zero
  ! them before the real transform.
  subroutine zero_vvdx_padding(m)
    integer(C_INT), intent(in) :: m
    integer(C_INT) :: i, j, k, y_first, y_last
    y_first = ny0 - 2
    y_last = nyN + 2
    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(VVdx, nx, nxd, nzB, y_first, y_last, m) private(i, j, k)
    do i = y_first, y_last
      do j = 1, nzB
        do k = nx + 2, nxd + 1
          VVdx(k, j, i, m) = 0.0d0
        end do
      end do
    end do
  end subroutine zero_vvdx_padding

  subroutine transform_to_physical()
    integer(C_INT) :: m
    call assemble_vvdz()
    do m = 1, 3
      call IFT(m)
      call transpose_zTOx_start(VVdz(:, :, :, m), m)
      if (m > 1) call to_x_points(m - 1)
    end do
    call to_x_points(3)
  end subroutine transform_to_physical

  subroutine to_x_points(m)
    integer(C_INT), intent(in) :: m
    call transpose_zTOx_finish(VVdz(:, :, :, m), VVdx(:, :, :, m), m)
    call zero_vvdx_padding(m)
    call RFT(m)
  end subroutine to_x_points

  ! products(:, :, :, 1:3) = uu, vv, ww (group 1) or uv, vw, uw (group 2),
  ! times factor = 1/(2 nxd nzd), the normalisation of the inverse
  ! transforms.  One pass over the three velocity fields: each thread reads
  ! u, v, w once and writes its three products (with the product index as
  ! the outer loop each field was streamed twice per call).
  subroutine build_products(g)
    integer(C_INT), intent(in) :: g
    integer(C_INT) :: i, j, k, y_first, y_last
    real(C_DOUBLE) :: u, v, w
    y_first = ny0 - 2
    y_last = nyN + 2
    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(rVVdx, products, nxd, nzB, y_first, y_last, factor, g) private(i, j, k, u, v, w)
    do i = y_first, y_last
      do j = 1, nzB
        do k = 1, 2*nxd
          u = rVVdx(k, j, i, 1); v = rVVdx(k, j, i, 2); w = rVVdx(k, j, i, 3)
          if (g == 1) then
            products(k, j, i, 1) = u*u*factor
            products(k, j, i, 2) = v*v*factor
            products(k, j, i, 3) = w*w*factor
          else
            products(k, j, i, 1) = u*v*factor
            products(k, j, i, 2) = v*w*factor
            products(k, j, i, 3) = u*w*factor
          end if
        end do
      end do
    end do
  end subroutine build_products

  ! products (physical) -> VVdz (spectral z-pencil).  The result for mode
  ! (iz, ix) of product p is VVdz(izd(iz) + 1, ix - nx0 + 1, iy, p).
  subroutine products_to_spectral()
    integer(C_INT) :: m
    do m = 1, 3
      call HFT(m)
      call transpose_xTOz_start(VVdp(:, :, :, m), m)
      if (m > 1) call to_z_modes(m - 1)
    end do
    call to_z_modes(3)
  end subroutine products_to_spectral

  subroutine to_z_modes(m)
    integer(C_INT), intent(in) :: m
    call transpose_xTOz_finish(VVdp(:, :, :, m), VVdz(:, :, :, m), m)
    call FFT(m)
  end subroutine to_z_modes

  ! Unpack VVdz(:, :, :, m) into a field with the layout of V(:, :, :, 1).
  subroutine vvdz_to_field(field, m)
    complex(C_DOUBLE_COMPLEX), intent(out) :: field(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
    integer(C_INT), intent(in) :: m
    integer(C_INT) :: ix, iy, iz, y_first, y_last
    y_first = ny0 - 2
    y_last = nyN + 2
    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(VVdz, izd, field, nx0, nxN, nz, y_first, y_last, m) private(ix, iy, iz)
    do ix = nx0, nxN
      do iy = y_first, y_last
        do iz = -nz, nz
          field(iy, iz, ix) = VVdz(izd(iz) + 1, ix - nx0 + 1, iy, m)
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
