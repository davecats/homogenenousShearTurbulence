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

    allocate (y(-2:ny + 1))
    do iy = -2, ny + 1
      y(iy) = iy*ly/ny
    end do
    allocate (izd(-nz:nz), ialfa(nx0:nxN), ibeta(-nz:nz), k2(-nz:nz, nx0:nxN))
    izd = [(merge(iz, nzd + iz, iz >= 0), iz=-nz, nz)]
    ialfa = [(dcmplx(0.0d0, ix*alfa0), ix=nx0, nxN)]
    ibeta = [(dcmplx(0.0d0, iz*beta0), iz=-nz, nz)]
    do ix = nx0, nxN
      do iz = -nz, nz
        k2(iz, ix) = (alfa0*ix)**2 + (beta0*iz)**2
      end do
    end do
    !$omp target enter data map(to: y, izd, ialfa, ibeta, k2, RK_rai)

    allocate (V(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN, 1:3)); V = 0
    allocate (oldrhs(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN, 1:2)); oldrhs = 0
    allocate (memrhs(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN, 1:2)); memrhs = 0
    !$omp target enter data map(to: V, oldrhs, memrhs)
  end subroutine allocate_fields

  subroutine free_fields()
    !$omp target exit data map(delete: V, oldrhs, memrhs, y, izd, ialfa, ibeta, k2, RK_rai)
    deallocate (V, oldrhs, memrhs, y, izd, ialfa, ibeta, k2)
  end subroutine free_fields

end module hst_setup
