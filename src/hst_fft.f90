! Fourier transforms in z (complex, in place) and in x (real <-> complex),
! batched over all y rows, and the four buffers they run in.
!
! Two backends, chosen at compile time:
!   HAVE_FFTW  FFTW3 on the host
!   HAVE_CUDA  cuFFT on the device; the buffers live on the device and the
!              plans run on their device addresses
! This file and hst_mpi.f90 are the only ones that name a vendor library.
! An AMD backend (hipFFT) goes here as a third block.
!
!   VVdz(iz, ix, iy)      z-pencil, padded to nzd, complex
!   VVdx(ix, iz, iy)      x-pencil, nxd+1 complex x modes
!   rVVdx(ix, iz, iy, 3)  physical space, 2*(nxd+1) reals in x (padded), u v w
!   products(ix, iz, iy)  one product of two velocity components, physical
!
! From channel/src/fft/ffts.fypp with HIP, the byte workspace and the
! double buffers for overlapped communication removed.
module hst_fft

  use, intrinsic :: iso_c_binding
  use hst_params
#ifdef HAVE_CUDA
  use cudafor
  use cufft
#endif

  implicit none
  private

  public :: init_fft, free_fft, FFT, IFT, RFT, HFT
  public :: VVdz, VVdx, rVVdx, products

  complex(C_DOUBLE_COMPLEX), allocatable, target, save :: VVdz(:, :, :), VVdx(:, :, :)
  real(C_DOUBLE), allocatable, target, save :: rVVdx(:, :, :, :), products(:, :, :)
  integer(C_INT), save :: fft_y0, fft_yN, fft_ny

#ifdef HAVE_FFTW
  include 'fftw3.f03'
  integer, save :: plan_type = FFTW_PATIENT
  type(C_PTR), save :: pFFT, pIFT, pRFT, pHFT
#endif
#ifdef HAVE_CUDA
  integer, save :: cu_pFFT, cu_pIFT, cu_pRFT, cu_pHFT
#endif

contains

  subroutine init_fft()
    integer :: istat
#ifdef HAVE_FFTW
    integer(C_INT), dimension(1) :: n_z, n_x, rn_x
#endif
#ifdef HAVE_CUDA
    integer, dimension(1) :: n, inembed, onembed
#endif
    fft_y0 = ny0 - 2
    fft_yN = nyN + 2
    fft_ny = fft_yN - fft_y0 + 1

    allocate (VVdz(nzd, nxB, fft_y0:fft_yN))
    allocate (VVdx(nxd + 1, nzB, fft_y0:fft_yN))
    allocate (rVVdx(2*(nxd + 1), nzB, fft_y0:fft_yN, 3))
    allocate (products(2*(nxd + 1), nzB, fft_y0:fft_yN))
    VVdz = 0; VVdx = 0; rVVdx = 0; products = 0
    !$omp target enter data map(to: VVdz, VVdx, rVVdx, products)

#ifdef HAVE_FFTW
    n_z = [nzd]; n_x = [nxd]; rn_x = [2*nxd]
    pFFT = fftw_plan_many_dft(1, n_z, nxB, VVdz(:, :, fft_y0), n_z, 1, nzd, &
                              VVdz(:, :, fft_y0), n_z, 1, nzd, FFTW_FORWARD, plan_type)
    pIFT = fftw_plan_many_dft(1, n_z, nxB, VVdz(:, :, fft_y0), n_z, 1, nzd, &
                              VVdz(:, :, fft_y0), n_z, 1, nzd, FFTW_BACKWARD, plan_type)
    pRFT = fftw_plan_many_dft_c2r(1, rn_x, nzB, VVdx(:, :, fft_y0), n_x + 1, 1, nxd + 1, &
                                  rVVdx(:, :, fft_y0, 1), 2*(n_x + 1), 1, 2*(nxd + 1), plan_type)
    pHFT = fftw_plan_many_dft_r2c(1, rn_x, nzB, rVVdx(:, :, fft_y0, 1), 2*(n_x + 1), 1, 2*(nxd + 1), &
                                  VVdx(:, :, fft_y0), n_x + 1, 1, nxd + 1, plan_type)
    istat = 0
#endif
#ifdef HAVE_CUDA
    istat = cufftPlan1d(cu_pIFT, nzd, CUFFT_Z2Z, fft_ny*nxB)
    call check(istat, 'cufftPlan1d IFT')
    istat = cufftPlan1d(cu_pFFT, nzd, CUFFT_Z2Z, fft_ny*nxB)
    call check(istat, 'cufftPlan1d FFT')
    n(1) = 2*nxd
    inembed(1) = nxd + 1
    onembed(1) = 2*(nxd + 1)
    istat = cufftPlanMany(cu_pRFT, 1, n, inembed, 1, nxd + 1, onembed, 1, 2*(nxd + 1), CUFFT_Z2D, nzB*fft_ny)
    call check(istat, 'cufftPlanMany RFT')
    istat = cufftPlanMany(cu_pHFT, 1, n, onembed, 1, 2*(nxd + 1), inembed, 1, nxd + 1, CUFFT_D2Z, nzB*fft_ny)
    call check(istat, 'cufftPlanMany HFT')
#endif
  end subroutine init_fft

  subroutine free_fft()
    integer :: istat
#ifdef HAVE_FFTW
    call fftw_destroy_plan(pFFT); call fftw_destroy_plan(pIFT)
    call fftw_destroy_plan(pRFT); call fftw_destroy_plan(pHFT)
    istat = 0
#endif
#ifdef HAVE_CUDA
    istat = cufftDestroy(cu_pFFT); istat = cufftDestroy(cu_pIFT)
    istat = cufftDestroy(cu_pRFT); istat = cufftDestroy(cu_pHFT)
#endif
    !$omp target exit data map(delete: VVdz, VVdx, rVVdx, products)
    deallocate (VVdz, VVdx, rVVdx, products)
  end subroutine free_fft

  subroutine check(istat, where)
    integer, intent(in) :: istat
    character(len=*), intent(in) :: where
    if (istat /= 0) then
      print *, trim(where), ' failed:', istat
      error stop 1
    end if
  end subroutine check

  ! Complex transform along z, in place, forward (FFT) or backward (IFT).
  subroutine FFT(x)
    complex(C_DOUBLE_COMPLEX), intent(inout), target :: x(:, :, ny0 - 2:)
    integer :: i, istat, y0
    y0 = lbound(x, 3)
#ifdef HAVE_FFTW
    do i = fft_y0, fft_yN
      call fftw_execute_dft(pFFT, x(:, :, i), x(:, :, i))
    end do
#endif
#ifdef HAVE_CUDA
    !$omp target data use_device_addr(x)
    istat = cudaDeviceSynchronize()
    istat = cufftExecZ2Z(cu_pFFT, x(1, 1, y0), x(1, 1, y0), CUFFT_FORWARD)
    call check(istat, 'cufftExecZ2Z FFT')
    istat = cudaDeviceSynchronize()
    !$omp end target data
#endif
  end subroutine FFT

  subroutine IFT(x)
    complex(C_DOUBLE_COMPLEX), intent(inout), target :: x(:, :, ny0 - 2:)
    integer :: i, istat, y0
    y0 = lbound(x, 3)
#ifdef HAVE_FFTW
    do i = fft_y0, fft_yN
      call fftw_execute_dft(pIFT, x(:, :, i), x(:, :, i))
    end do
#endif
#ifdef HAVE_CUDA
    !$omp target data use_device_addr(x)
    istat = cudaDeviceSynchronize()
    istat = cufftExecZ2Z(cu_pIFT, x(1, 1, y0), x(1, 1, y0), CUFFT_INVERSE)
    call check(istat, 'cufftExecZ2Z IFT')
    istat = cudaDeviceSynchronize()
    !$omp end target data
#endif
  end subroutine IFT

  ! Complex x modes -> real x points (RFT) and back (HFT).
  subroutine RFT(x, rx)
    complex(C_DOUBLE_COMPLEX), intent(inout), target :: x(:, :, ny0 - 2:)
    real(C_DOUBLE), intent(inout), target :: rx(:, :, ny0 - 2:)
    integer :: i, istat, x_y0, rx_y0
    x_y0 = lbound(x, 3)
    rx_y0 = lbound(rx, 3)
#ifdef HAVE_FFTW
    do i = fft_y0, fft_yN
      call fftw_execute_dft_c2r(pRFT, x(:, :, i), rx(:, :, i))
    end do
#endif
#ifdef HAVE_CUDA
    !$omp target data use_device_addr(x, rx)
    istat = cudaDeviceSynchronize()
    istat = cufftExecZ2D(cu_pRFT, x(1, 1, x_y0), rx(1, 1, rx_y0))
    call check(istat, 'cufftExecZ2D RFT')
    istat = cudaDeviceSynchronize()
    !$omp end target data
#endif
  end subroutine RFT

  subroutine HFT(rx, x)
    real(C_DOUBLE), intent(inout), target :: rx(:, :, ny0 - 2:)
    complex(C_DOUBLE_COMPLEX), intent(inout), target :: x(:, :, ny0 - 2:)
    integer :: i, istat, x_y0, rx_y0
    x_y0 = lbound(x, 3)
    rx_y0 = lbound(rx, 3)
#ifdef HAVE_FFTW
    do i = fft_y0, fft_yN
      call fftw_execute_dft_r2c(pHFT, rx(:, :, i), x(:, :, i))
    end do
#endif
#ifdef HAVE_CUDA
    !$omp target data use_device_addr(rx, x)
    istat = cudaDeviceSynchronize()
    istat = cufftExecD2Z(cu_pHFT, rx(1, 1, rx_y0), x(1, 1, x_y0))
    call check(istat, 'cufftExecD2Z HFT')
    istat = cudaDeviceSynchronize()
    !$omp end target data
#endif
  end subroutine HFT

end module hst_fft
