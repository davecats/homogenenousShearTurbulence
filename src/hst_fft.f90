! Fourier transforms in z (complex, in place) and in x (real <-> complex, in
! place), batched over all y rows and over the three fields of a buffer, and
! the buffers they run in.
!
! Two backends, chosen at compile time:
!   HAVE_FFTW  FFTW3 on the host
!   HAVE_CUDA  cuFFT on the device; the buffers live on the device and the
!              plans run on their device addresses
! This file and hst_mpi.f90 are the only ones that name a vendor library.
! An AMD backend (hipFFT) goes here as a third block.
!
!   VVdz(iz, ix, iy, m)     z-pencil, padded to nzd, complex; m = u, v, w or three products
!   VVdx(ix, iz, iy, m)     x-pencil, nxd+1 complex x modes of u, v, w ...
!   rVVdx(ix, iz, iy, m)    ... and the same bytes seen as 2(nxd+1) real x points
!   VVdp, products          the same pair for three products
!
! The x transforms run in place: a row of nxd+1 complex modes and a row of
! 2(nxd+1) reals are the same bytes (the padded in-place layout of FFTW and
! cuFFT), and the real names are pointer views of the complex buffers.  On
! the device a view is found through the mapping of the buffer it points
! into.  Three fields per transform and per transpose: the u, v, w of a
! substep, or three of the six products.
!
! From channel/src/fft/ffts.fypp with HIP, the byte workspace and the
! double buffers for overlapped communication removed.
module hst_fft

  use, intrinsic :: iso_c_binding
  use hst_params
#ifdef HAVE_CUDA
  use cudafor
  use cufft
  use omp_lib, only: omp_get_default_device
#endif

  implicit none
  private

  public :: init_fft, free_fft, FFT, IFT, RFT, HFT, device_sync
  public :: VVdz, VVdx, rVVdx, VVdp, products

  complex(C_DOUBLE_COMPLEX), allocatable, target, save :: VVdz(:, :, :, :), VVdx(:, :, :, :), VVdp(:, :, :, :)
  real(C_DOUBLE), pointer, save :: rVVdx(:, :, :, :), products(:, :, :, :)
  integer(C_INT), save :: fft_y0, fft_yN, fft_ny

#ifdef HAVE_FFTW
  include 'fftw3.f03'
  integer, save :: plan_type = FFTW_MEASURE      ! PATIENT plans take minutes at 256^3 for a few % of transform time
  type(C_PTR), save :: pFFT, pIFT, pRFT, pHFT
#endif
#ifdef HAVE_CUDA
  integer, save :: cu_pFFT, cu_pIFT, cu_pRFT, cu_pHFT
  ! One work area for the four plans (they never run at the same time):
  ! cuFFT would otherwise give each plan its own, about the size of the
  ! data it transforms, which at 512^3 on one GPU is more than the fields.
  integer(1), device, allocatable, save :: work(:)
  ! The CUDA stream the OpenMP target regions run on (NVHPC extension).
  ! The cuFFT plans are put on it, so transforms and kernels are ordered
  ! by the stream and no host synchronisation is needed between them.
  interface
    function ompx_get_cuda_stream(device, nowait) bind(c, name='ompx_get_cuda_stream') result(stream)
      import :: C_PTR, C_INT
      integer(C_INT), value :: device, nowait
      type(C_PTR) :: stream
    end function ompx_get_cuda_stream
  end interface
#endif

contains

  subroutine init_fft()
    integer :: istat
#ifdef HAVE_FFTW
    integer(C_INT), dimension(1) :: n_z, n_x, e_c, e_r
#endif
#ifdef HAVE_CUDA
    integer, dimension(1) :: n, inembed, onembed
    integer(kind=cuda_stream_kind) :: stream
    integer(C_SIZE_T) :: ws(4)
#endif
    fft_y0 = ny0 - 2
    fft_yN = nyN + 2
    fft_ny = fft_yN - fft_y0 + 1

    allocate (VVdz(nzd, nxB, fft_y0:fft_yN, 3))
    allocate (VVdx(nxd + 1, nzB, fft_y0:fft_yN, 3), VVdp(nxd + 1, nzB, fft_y0:fft_yN, 3))
    VVdz = 0; VVdx = 0; VVdp = 0
    call c_f_pointer(c_loc(VVdx), rVVdx, [2*(nxd + 1), nzB, fft_ny, 3])
    call c_f_pointer(c_loc(VVdp), products, [2*(nxd + 1), nzB, fft_ny, 3])
    rVVdx(1:, 1:, fft_y0:, 1:) => rVVdx
    products(1:, 1:, fft_y0:, 1:) => products
    !$omp target enter data map(to: VVdz, VVdx, VVdp)

#ifdef HAVE_FFTW
    n_z = [nzd]; n_x = [2*nxd]; e_c = [nxd + 1]; e_r = [2*(nxd + 1)]
    pFFT = fftw_plan_many_dft(1, n_z, nxB*fft_ny*3, VVdz, n_z, 1, nzd, VVdz, n_z, 1, nzd, FFTW_FORWARD, plan_type)
    pIFT = fftw_plan_many_dft(1, n_z, nxB*fft_ny*3, VVdz, n_z, 1, nzd, VVdz, n_z, 1, nzd, FFTW_BACKWARD, plan_type)
    pRFT = fftw_plan_many_dft_c2r(1, n_x, nzB*fft_ny*3, VVdx, e_c, 1, nxd + 1, rVVdx, e_r, 1, 2*(nxd + 1), plan_type)
    pHFT = fftw_plan_many_dft_r2c(1, n_x, nzB*fft_ny*3, products, e_r, 1, 2*(nxd + 1), VVdp, e_c, 1, nxd + 1, plan_type)
    istat = 0
#endif
#ifdef HAVE_CUDA
    n(1) = 2*nxd
    inembed(1) = nxd + 1
    onembed(1) = 2*(nxd + 1)
    istat = cufftCreate(cu_pIFT); istat = cufftSetAutoAllocation(cu_pIFT, 0)
    istat = cufftMakePlan1d(cu_pIFT, nzd, CUFFT_Z2Z, fft_ny*nxB*3, ws(1))
    call check(istat, 'cufftMakePlan1d IFT')
    istat = cufftCreate(cu_pFFT); istat = cufftSetAutoAllocation(cu_pFFT, 0)
    istat = cufftMakePlan1d(cu_pFFT, nzd, CUFFT_Z2Z, fft_ny*nxB*3, ws(2))
    call check(istat, 'cufftMakePlan1d FFT')
    istat = cufftCreate(cu_pRFT); istat = cufftSetAutoAllocation(cu_pRFT, 0)
    istat = cufftMakePlanMany(cu_pRFT, 1, n, inembed, 1, nxd + 1, onembed, 1, 2*(nxd + 1), CUFFT_Z2D, nzB*fft_ny*3, ws(3))
    call check(istat, 'cufftMakePlanMany RFT')
    istat = cufftCreate(cu_pHFT); istat = cufftSetAutoAllocation(cu_pHFT, 0)
    istat = cufftMakePlanMany(cu_pHFT, 1, n, onembed, 1, 2*(nxd + 1), inembed, 1, nxd + 1, CUFFT_D2Z, nzB*fft_ny*3, ws(4))
    call check(istat, 'cufftMakePlanMany HFT')
    allocate (work(max(maxval(ws), 1_C_SIZE_T)))
    istat = cufftSetWorkArea(cu_pIFT, work); call check(istat, 'cufftSetWorkArea IFT')
    istat = cufftSetWorkArea(cu_pFFT, work); call check(istat, 'cufftSetWorkArea FFT')
    istat = cufftSetWorkArea(cu_pRFT, work); call check(istat, 'cufftSetWorkArea RFT')
    istat = cufftSetWorkArea(cu_pHFT, work); call check(istat, 'cufftSetWorkArea HFT')
    if (has_terminal) write (*, '(A,F8.1,A)') '   cuFFT work area: ', maxval(ws)/1024.0d0**2, ' MB per rank'
    stream = transfer(ompx_get_cuda_stream(int(omp_get_default_device(), C_INT), 0_C_INT), stream)
    istat = cufftSetStream(cu_pFFT, stream); call check(istat, 'cufftSetStream FFT')
    istat = cufftSetStream(cu_pIFT, stream); call check(istat, 'cufftSetStream IFT')
    istat = cufftSetStream(cu_pRFT, stream); call check(istat, 'cufftSetStream RFT')
    istat = cufftSetStream(cu_pHFT, stream); call check(istat, 'cufftSetStream HFT')
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
    deallocate (work)
#endif
    !$omp target exit data map(delete: VVdz, VVdx, VVdp)
    nullify (rVVdx, products)
    deallocate (VVdz, VVdx, VVdp)
  end subroutine free_fft

  ! Wait for everything queued on the device (the timer's boundaries).
  subroutine device_sync()
#ifdef HAVE_CUDA
    integer :: istat
    istat = cudaDeviceSynchronize()
#endif
  end subroutine device_sync

  subroutine check(istat, where)
    integer, intent(in) :: istat
    character(len=*), intent(in) :: where
    if (istat /= 0) then
      print *, trim(where), ' failed:', istat
      error stop 1
    end if
  end subroutine check

  ! Complex transform of VVdz along z, in place: forward (FFT) or backward (IFT).
  subroutine FFT()
    integer :: istat
#ifdef HAVE_FFTW
    call fftw_execute_dft(pFFT, VVdz, VVdz)
#endif
#ifdef HAVE_CUDA
    !$omp target data use_device_addr(VVdz)
    istat = cufftExecZ2Z(cu_pFFT, VVdz, VVdz, CUFFT_FORWARD)
    call check(istat, 'cufftExecZ2Z FFT')
    !$omp end target data
#endif
  end subroutine FFT

  subroutine IFT()
    integer :: istat
#ifdef HAVE_FFTW
    call fftw_execute_dft(pIFT, VVdz, VVdz)
#endif
#ifdef HAVE_CUDA
    !$omp target data use_device_addr(VVdz)
    istat = cufftExecZ2Z(cu_pIFT, VVdz, VVdz, CUFFT_INVERSE)
    call check(istat, 'cufftExecZ2Z IFT')
    !$omp end target data
#endif
  end subroutine IFT

  ! Complex x modes of VVdx -> real x points rVVdx (RFT), and the real
  ! products -> complex modes VVdp (HFT), both in place.
  subroutine RFT()
    integer :: istat
#ifdef HAVE_FFTW
    call fftw_execute_dft_c2r(pRFT, VVdx, rVVdx)
#endif
#ifdef HAVE_CUDA
    !$omp target data use_device_addr(VVdx)
    istat = cufftExecZ2D(cu_pRFT, VVdx, VVdx)
    call check(istat, 'cufftExecZ2D RFT')
    !$omp end target data
#endif
  end subroutine RFT

  subroutine HFT()
    integer :: istat
#ifdef HAVE_FFTW
    call fftw_execute_dft_r2c(pHFT, products, VVdp)
#endif
#ifdef HAVE_CUDA
    !$omp target data use_device_addr(VVdp)
    istat = cufftExecD2Z(cu_pHFT, VVdp, VVdp)
    call check(istat, 'cufftExecD2Z HFT')
    !$omp end target data
#endif
  end subroutine HFT

end module hst_fft
