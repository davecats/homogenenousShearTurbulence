! The x-z pencil decomposition and the transposes between its two layouts.
!
! In spectral space (z-pencil) a rank owns the x modes nx0:nxN and every z
! mode; in physical space (x-pencil) it owns the z lines nz0:nzN and every x
! point.  Moving between the two is one alltoall over all ranks, carrying
! three fields at once (u, v, w, or three products).  Each rank owns the
! whole of y (npy = 1); the ny0:nyN names are kept so that a y
! decomposition can be added later.
!
! The alltoall itself goes through MPI (CUDA-aware, the buffers stay on the
! device) or, in a build with NCCL=1 and one GPU per rank, through NCCL as
! grouped send/recv pairs on the OpenMP target stream, so that no host
! synchronisation is needed around it (deck parameter transport).
!
! Taken from channel/src/mpi/mpi_transpose.f90 with the y-slab machinery
! and HIP removed.  The pack/unpack kernels keep two rules from there:
! the innermost loop is the index the *read* runs contiguously in, and
! pack_zTOx/unpack_zTOx (likewise pack_xTOz/unpack_xTOz) spell the buffer
! position p identically, because the alltoall permutes whole blocks.
module hst_mpi

  use, intrinsic :: iso_c_binding
  use mpi_f08
  use hst_params
  use hst_timer, only: toc, T_PACK, T_ALLTOALL
#ifdef HAVE_NCCL
  use hst_fft, only: target_stream
  use omp_lib, only: omp_get_num_devices, omp_get_default_device
  use cudafor, only: cudaSetDevice
#endif

  implicit none
  private

  public :: setup_decomposition, free_mpi, transpose_zTOx, transpose_xTOz
  public :: cpl_view_type, cpl_pview_type

  complex(C_DOUBLE_COMPLEX), allocatable, target, save :: sendbuf(:), recvbuf(:)
  integer(C_INT), save :: sendcount
  !$omp declare target(sendcount)
  logical, save :: transpose_is_local
  logical, save :: use_nccl = .false.

#ifdef HAVE_NCCL
  ! NCCL through its C prototypes (nccl.h).  The Fortran module of NVHPC
  ! wants CUDA Fortran device arrays, which OpenMP-mapped arrays are not.
  type, bind(c) :: nccl_unique_id
    integer(C_INT8_T) :: bytes(128)
  end type nccl_unique_id
  integer(C_INT), parameter :: NCCL_UINT8 = 1
  type(C_PTR), save :: nccl_comm = C_NULL_PTR, nccl_stream = C_NULL_PTR
  interface
    function ncclGetUniqueId(id) bind(c, name='ncclGetUniqueId') result(r)
      import :: nccl_unique_id, C_INT
      type(nccl_unique_id) :: id
      integer(C_INT) :: r
    end function ncclGetUniqueId
    function ncclCommInitRank(comm, nranks, id, rank) bind(c, name='ncclCommInitRank') result(r)
      import :: nccl_unique_id, C_INT, C_PTR
      type(C_PTR) :: comm
      integer(C_INT), value :: nranks, rank
      type(nccl_unique_id), value :: id
      integer(C_INT) :: r
    end function ncclCommInitRank
    function ncclCommDestroy(comm) bind(c, name='ncclCommDestroy') result(r)
      import :: C_INT, C_PTR
      type(C_PTR), value :: comm
      integer(C_INT) :: r
    end function ncclCommDestroy
    function ncclGroupStart() bind(c, name='ncclGroupStart') result(r)
      import :: C_INT
      integer(C_INT) :: r
    end function ncclGroupStart
    function ncclGroupEnd() bind(c, name='ncclGroupEnd') result(r)
      import :: C_INT
      integer(C_INT) :: r
    end function ncclGroupEnd
    function ncclSend(buf, count, dtype, peer, comm, stream) bind(c, name='ncclSend') result(r)
      import :: C_INT, C_PTR, C_SIZE_T
      type(C_PTR), value :: buf, comm, stream
      integer(C_SIZE_T), value :: count
      integer(C_INT), value :: dtype, peer
      integer(C_INT) :: r
    end function ncclSend
    function ncclRecv(buf, count, dtype, peer, comm, stream) bind(c, name='ncclRecv') result(r)
      import :: C_INT, C_PTR, C_SIZE_T
      type(C_PTR), value :: buf, comm, stream
      integer(C_SIZE_T), value :: count
      integer(C_INT), value :: dtype, peer
      integer(C_INT) :: r
    end function ncclRecv
  end interface
#endif
  ! MPI-IO views of this rank's x slab in a CPL-layout file (see hst_io):
  ! the file array is (3, ny+4, 2nz+1, nx+1) in Fortran order for the
  ! velocity and (ny+4, 2nz+1, nx+1) for the pressure.
  type(MPI_Datatype), save :: cpl_view_type, cpl_pview_type
  integer :: ierr

contains

  ! npxz = nproc ranks each own nxB = (nx+1)/nproc x modes in spectral space
  ! and nzB = nzd/nproc z lines in physical space.  The alltoall needs both
  ! splits to be even.
  subroutine setup_decomposition()
    integer(C_SIZE_T) :: n

    npxz = nproc
    ipxz = iproc
    if (mod(nx + 1, npxz) /= 0 .or. mod(nzd, npxz) /= 0) then
      if (has_terminal) then
        print *, 'ERROR: nproc must divide both nx+1 and nzd.'
        print *, '       nx+1 =', nx + 1, ' nzd =', nzd, ' nproc =', nproc
      end if
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end if
    nx0 = ipxz*(nx + 1)/npxz
    nxN = (ipxz + 1)*(nx + 1)/npxz - 1
    nxB = nxN - nx0 + 1
    nz0 = ipxz*nzd/npxz
    nzN = (ipxz + 1)*nzd/npxz - 1
    nzB = nzN - nz0 + 1
    ny0 = 0
    nyN = ny - 1
    has_average = (nx0 == 0)
    !$omp target update to(nx0, nxN, nxB, nz0, nzN, nzB, ny0, nyN, ny, ni, S)
    if (has_terminal) write (*, '(A,I5,A,I5,A,I5)') '   ranks =', nproc, '   nxB   =', nxB, '   nzB   =', nzB

    transpose_is_local = (nzB == nzd)
    sendcount = nxB*nzB*(nyN - ny0 + 5)*3          ! three fields per transpose
    !$omp target update to(sendcount)
    n = 1
    if (.not. transpose_is_local) n = int(npxz, C_SIZE_T)*int(sendcount, C_SIZE_T)
    allocate (sendbuf(n), recvbuf(n))
    sendbuf = 0; recvbuf = 0
    !$omp target enter data map(alloc: sendbuf, recvbuf)
    call setup_transport()

    call MPI_Type_create_subarray(4, [3, ny + 4, 2*nz + 1, nx + 1], [3, ny + 4, 2*nz + 1, nxB], &
                                  [0, 0, 0, nx0], MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, cpl_view_type, ierr)
    call MPI_Type_commit(cpl_view_type, ierr)
    call MPI_Type_create_subarray(3, [ny + 4, 2*nz + 1, nx + 1], [ny + 4, 2*nz + 1, nxB], &
                                  [0, 0, nx0], MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, cpl_pview_type, ierr)
    call MPI_Type_commit(cpl_pview_type, ierr)
  end subroutine setup_decomposition

  subroutine free_mpi()
#ifdef HAVE_NCCL
    if (use_nccl) ierr = ncclCommDestroy(nccl_comm)
#endif
    !$omp target exit data map(delete: sendbuf, recvbuf)
    deallocate (sendbuf, recvbuf)
    call MPI_Type_free(cpl_view_type, ierr)
    call MPI_Type_free(cpl_pview_type, ierr)
  end subroutine free_mpi

  ! transport = 'nccl' needs a build with NCCL=1 and one GPU per rank;
  ! 'auto' takes NCCL when both hold and MPI otherwise.  One rank needs no
  ! transport at all.
  subroutine setup_transport()
#ifdef HAVE_NCCL
    type(nccl_unique_id) :: id
    type(MPI_Comm) :: node
    integer :: node_ranks, r
#endif
    use_nccl = .false.
    if (transpose_is_local) return
    if (transport /= 'mpi') then
#ifdef HAVE_NCCL
      call MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, node, ierr)
      call MPI_Comm_size(node, node_ranks, ierr)
      call MPI_Comm_free(node, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, node_ranks, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
      if (node_ranks <= omp_get_num_devices()) then
        if (iproc == 0) r = ncclGetUniqueId(id)
        call MPI_Bcast(id%bytes, 128, MPI_BYTE, 0, MPI_COMM_WORLD, ierr)
        r = cudaSetDevice(omp_get_default_device())
        r = ncclCommInitRank(nccl_comm, nproc, id, iproc)
        if (r /= 0) then
          if (has_terminal) print *, 'ERROR: ncclCommInitRank failed with code', r
          call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        nccl_stream = target_stream()
        use_nccl = .true.
      else if (transport == 'nccl') then
        if (has_terminal) print *, 'ERROR: transport = nccl needs one GPU per rank'
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
      end if
#else
      if (transport == 'nccl') then
        if (has_terminal) print *, 'ERROR: transport = nccl needs a build with NCCL=1'
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
      end if
#endif
    end if
    if (has_terminal) write (*, '(A,A)') '   alltoall transport: ', merge('nccl', 'mpi ', use_nccl)
  end subroutine setup_transport

  !------------------------------------------------------------------------
  ! z-pencil Vz(iz, ix, iy, m)  <->  x-pencil Vx(ix, iz, iy, m), m = 1..3
  !------------------------------------------------------------------------

  subroutine repack_zTOx_local(Vz, Vx)
    complex(C_DOUBLE_COMPLEX), intent(in) :: Vz(:, :, :, :)
    complex(C_DOUBLE_COMPLEX), intent(out) :: Vx(:, :, :, :)
    integer(C_SIZE_T) :: iy, ix, iz, m
    integer(C_INT) :: ny_batch, nb
    ny_batch = size(Vz, 3); nb = size(Vz, 4)
    !$omp target teams distribute parallel do collapse(4) default(none) &
    !$omp shared(Vz, Vx, ny_batch, nb, nxB, nzd) private(iy, ix, iz, m)
    do m = 1, nb
      do iy = 1, ny_batch
        do ix = 1, nxB
          do iz = 1, nzd
            Vx(ix, iz, iy, m) = Vz(iz, ix, iy, m)
          end do
        end do
      end do
    end do
  end subroutine repack_zTOx_local

  subroutine repack_xTOz_local(Vx, Vz)
    complex(C_DOUBLE_COMPLEX), intent(in) :: Vx(:, :, :, :)
    complex(C_DOUBLE_COMPLEX), intent(out) :: Vz(:, :, :, :)
    integer(C_SIZE_T) :: iy, ix, iz, m
    integer(C_INT) :: ny_batch, nb
    ny_batch = size(Vx, 3); nb = size(Vx, 4)
    !$omp target teams distribute parallel do collapse(4) default(none) &
    !$omp shared(Vx, Vz, ny_batch, nb, nxB, nzd) private(iy, ix, iz, m)
    do m = 1, nb
      do iy = 1, ny_batch
        do iz = 1, nzd
          do ix = 1, nxB
            Vz(iz, ix, iy, m) = Vx(ix, iz, iy, m)
          end do
        end do
      end do
    end do
  end subroutine repack_xTOz_local

  subroutine pack_zTOx(Vz, send)
    complex(C_DOUBLE_COMPLEX), intent(in) :: Vz(:, :, :, :)
    complex(C_DOUBLE_COMPLEX), intent(out) :: send(:)
    integer(C_SIZE_T) :: iy, ix, iz, m, dest, p
    integer(C_INT) :: ny_batch, nb
    ny_batch = size(Vz, 3); nb = size(Vz, 4)
    !$omp target teams distribute parallel do collapse(5) default(none) &
    !$omp shared(Vz, send, ny_batch, nb, nxB, nzB, npxz, sendcount) private(iy, ix, iz, m, dest, p)
    do dest = 0, npxz - 1
      do m = 1, nb
        do iy = 1, ny_batch
          do ix = 1, nxB
            do iz = 1, nzB
              p = dest*sendcount + iz + nzB*(ix - 1) + nzB*nxB*(iy - 1) + nzB*nxB*ny_batch*(m - 1)
              send(p) = Vz(dest*nzB + iz, ix, iy, m)
            end do
          end do
        end do
      end do
    end do
  end subroutine pack_zTOx

  subroutine unpack_zTOx(recv, Vx)
    complex(C_DOUBLE_COMPLEX), intent(in) :: recv(:)
    complex(C_DOUBLE_COMPLEX), intent(out) :: Vx(:, :, :, :)
    integer(C_SIZE_T) :: iy, ix, iz, m, src, p
    integer(C_INT) :: ny_batch, nb
    ny_batch = size(Vx, 3); nb = size(Vx, 4)
    !$omp target teams distribute parallel do collapse(5) default(none) &
    !$omp shared(Vx, recv, ny_batch, nb, nxB, nzB, npxz, sendcount) private(iy, ix, iz, m, src, p)
    do src = 0, npxz - 1
      do m = 1, nb
        do iy = 1, ny_batch
          do ix = 1, nxB
            do iz = 1, nzB
              p = src*sendcount + iz + nzB*(ix - 1) + nzB*nxB*(iy - 1) + nzB*nxB*ny_batch*(m - 1)
              Vx(ix + src*nxB, iz, iy, m) = recv(p)
            end do
          end do
        end do
      end do
    end do
  end subroutine unpack_zTOx

  subroutine pack_xTOz(Vx, send)
    complex(C_DOUBLE_COMPLEX), intent(in) :: Vx(:, :, :, :)
    complex(C_DOUBLE_COMPLEX), intent(out) :: send(:)
    integer(C_SIZE_T) :: iy, ix, iz, m, dest, p
    integer(C_INT) :: ny_batch, nb
    ny_batch = size(Vx, 3); nb = size(Vx, 4)
    !$omp target teams distribute parallel do collapse(5) default(none) &
    !$omp shared(Vx, send, ny_batch, nb, nxB, nzB, npxz, sendcount) private(iy, ix, iz, m, dest, p)
    do dest = 0, npxz - 1
      do m = 1, nb
        do iy = 1, ny_batch
          do iz = 1, nzB
            do ix = 1, nxB
              p = dest*sendcount + ix + nxB*(iz - 1) + nxB*nzB*(iy - 1) + nxB*nzB*ny_batch*(m - 1)
              send(p) = Vx(dest*nxB + ix, iz, iy, m)
            end do
          end do
        end do
      end do
    end do
  end subroutine pack_xTOz

  subroutine unpack_xTOz(recv, Vz)
    complex(C_DOUBLE_COMPLEX), intent(in) :: recv(:)
    complex(C_DOUBLE_COMPLEX), intent(out) :: Vz(:, :, :, :)
    integer(C_SIZE_T) :: iy, ix, iz, m, src, p
    integer(C_INT) :: ny_batch, nb
    ny_batch = size(Vz, 3); nb = size(Vz, 4)
    !$omp target teams distribute parallel do collapse(5) default(none) &
    !$omp shared(Vz, recv, ny_batch, nb, nxB, nzB, npxz, sendcount) private(iy, ix, iz, m, src, p)
    do src = 0, npxz - 1
      do m = 1, nb
        do iy = 1, ny_batch
          do iz = 1, nzB
            do ix = 1, nxB
              p = src*sendcount + ix + nxB*(iz - 1) + nxB*nzB*(iy - 1) + nxB*nzB*ny_batch*(m - 1)
              Vz(iz + src*nzB, ix, iy, m) = recv(p)
            end do
          end do
        end do
      end do
    end do
  end subroutine unpack_xTOz

  ! The one collective of the solver.  On the GPU the buffers stay on the
  ! device: the use_device_addr block hands MPI (CUDA-aware) or NCCL their
  ! device addresses.  The NCCL form is one send and one receive per peer
  ! in a group (NCCL has no alltoall) on the OpenMP stream, so the host
  ! does not wait for it: the unpack kernel is queued behind it.
  subroutine alltoall()
#ifdef HAVE_NCCL
    integer(C_SIZE_T) :: nbytes
    integer :: peer, r
#endif
#ifdef HAVE_CUDA
    !$omp target data use_device_addr(sendbuf, recvbuf)
#endif
    if (use_nccl) then
#ifdef HAVE_NCCL
      nbytes = 16_C_SIZE_T*int(sendcount, C_SIZE_T)
      r = ncclGroupStart()
      do peer = 0, npxz - 1
        r = ncclSend(c_loc(sendbuf(peer*sendcount + 1)), nbytes, NCCL_UINT8, peer, nccl_comm, nccl_stream)
        if (r /= 0) error stop 'ncclSend failed'
        r = ncclRecv(c_loc(recvbuf(peer*sendcount + 1)), nbytes, NCCL_UINT8, peer, nccl_comm, nccl_stream)
        if (r /= 0) error stop 'ncclRecv failed'
      end do
      r = ncclGroupEnd()
      if (r /= 0) error stop 'ncclGroupEnd failed'
#endif
    else
      call MPI_Alltoall(sendbuf, int(sendcount), MPI_DOUBLE_COMPLEX, &
                        recvbuf, int(sendcount), MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)
      if (ierr /= MPI_SUCCESS) error stop 'MPI_Alltoall failed'
    end if
#ifdef HAVE_CUDA
    !$omp end target data
#endif
  end subroutine alltoall

  subroutine transpose_zTOx(Vz, Vx)
    complex(C_DOUBLE_COMPLEX), intent(in) :: Vz(:, :, :, :)
    complex(C_DOUBLE_COMPLEX), intent(out) :: Vx(:, :, :, :)
    if (transpose_is_local) then
      call repack_zTOx_local(Vz, Vx)
      call toc(T_PACK)
    else
      call pack_zTOx(Vz, sendbuf)
      call toc(T_PACK)
      call alltoall()
      call toc(T_ALLTOALL)
      call unpack_zTOx(recvbuf, Vx)
      call toc(T_PACK)
    end if
  end subroutine transpose_zTOx

  subroutine transpose_xTOz(Vx, Vz)
    complex(C_DOUBLE_COMPLEX), intent(in) :: Vx(:, :, :, :)
    complex(C_DOUBLE_COMPLEX), intent(out) :: Vz(:, :, :, :)
    if (transpose_is_local) then
      call repack_xTOz_local(Vx, Vz)
      call toc(T_PACK)
    else
      call pack_xTOz(Vx, sendbuf)
      call toc(T_PACK)
      call alltoall()
      call toc(T_ALLTOALL)
      call unpack_xTOz(recvbuf, Vz)
      call toc(T_PACK)
    end if
  end subroutine transpose_xTOz

end module hst_mpi
