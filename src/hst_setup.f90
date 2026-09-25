! Allocating the grid arrays and the fields, and mapping them to the device.
module hst_setup

  use, intrinsic :: iso_c_binding
  use hst_params

  implicit none
  private
  public :: allocate_fields, free_fields

contains

  subroutine allocate_fields()
    integer(C_INT) :: ix, iy, iz

    ! y grid: uniform, or clustered at mid-box with the tanh map of
    ! scddnsdata.cpl (htcoeff); the ghost rows are its periodic extension
    allocate (y(-2:ny + 1), dyl(0:ny - 1), fy(-2:ny + 1), inlayer(0:ny - 1))
    do iy = -2, ny + 1
      y(iy) = ly*ymap(real(modulo(iy, ny), C_DOUBLE)/ny) + ly*((iy - modulo(iy, ny))/ny)
    end do
    do iy = 0, ny - 1
      dyl(iy) = 0.5d0*(y(iy + 1) - y(iy - 1))
    end do
    fy = 0.0d0
    inlayer = 0
    allocate (izd(-nz:nz), ialfa(nx0:nxN), ibeta(-nz:nz), k2(-nz:nz, nx0:nxN))
    izd = [(merge(iz, nzd + iz, iz >= 0), iz=-nz, nz)]
    ialfa = [(dcmplx(0.0d0, ix*alfa0), ix=nx0, nxN)]
    ibeta = [(dcmplx(0.0d0, iz*beta0), iz=-nz, nz)]
    do ix = nx0, nxN
      do iz = -nz, nz
        k2(iz, ix) = (alfa0*ix)**2 + (beta0*iz)**2
      end do
    end do
    !$omp target enter data map(to: y, dyl, fy, inlayer, izd, ialfa, ibeta, k2, RK_rai)

    allocate (V(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN, 1:3)); V = 0
    allocate (oldrhs(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN, 1:2)); oldrhs = 0
    allocate (memrhs(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN, 1:2)); memrhs = 0
    !$omp target enter data map(to: V, oldrhs, memrhs)
  end subroutine allocate_fields

  subroutine free_fields()
    !$omp target exit data map(delete: V, oldrhs, memrhs, y, dyl, fy, inlayer, izd, ialfa, ibeta, k2, RK_rai)
    deallocate (V, oldrhs, memrhs, y, dyl, fy, inlayer, izd, ialfa, ibeta, k2)
  end subroutine free_fields

  ! xi in [0, 1) -> [0, 1): uniform for ystretch = 0, else the two-sided
  ! tanh clustering at 1/2 of scddnsdata.cpl.
  real(C_DOUBLE) function ymap(xi)
    real(C_DOUBLE), intent(in) :: xi
    if (ystretch <= 1.0d-10) then
      ymap = xi
    else if (xi <= 0.5d0) then
      ymap = 0.5d0*tanh(2.0d0*ystretch*xi)/tanh(ystretch)
    else
      ymap = 0.5d0*(2.0d0 + tanh(2.0d0*ystretch*(xi - 1.0d0))/tanh(ystretch))
    end if
  end function ymap

end module hst_setup
