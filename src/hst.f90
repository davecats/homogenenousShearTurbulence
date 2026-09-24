!============================================!
!                                            !
!     Direct Numerical Simulation (DNS)      !
!     of homogeneous shear turbulence        !
!                                            !
!============================================!
!
! Derived from the `channel` code (D. Gatti) and the CPL HST code.
! Usage:  mpirun -np N ./hst [hst.in]
!
program hst

  use, intrinsic :: iso_c_binding
  use mpi_f08
  use hst_params
  use hst_input
#ifdef HAVE_CUDA
  use omp_lib
#endif

  implicit none

  character(len=256) :: deck
  integer :: ierr

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, iproc, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
  has_terminal = (iproc == 0)

#ifdef HAVE_CUDA
  call select_device()
#endif

  deck = 'hst.in'
  if (command_argument_count() >= 1) call get_command_argument(1, deck)
  call read_input(trim(deck))
  if (has_terminal) call print_input()
  call setup_decomposition()

  ! WP1 and later: allocate fields, FFT plans, initial field, time loop.

  call MPI_Finalize(ierr)

contains

  ! 1-D x-z pencil decomposition: npxz = nproc ranks each own nxB = (nx+1)/nproc
  ! x modes in spectral space and nzB = nzd/nproc z lines in physical space,
  ! and the whole of y.  The alltoall transpose needs both splits to be even.
  subroutine setup_decomposition()
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
  end subroutine setup_decomposition

#ifdef HAVE_CUDA
  ! One GPU per rank: node-local rank modulo the number of devices.  The
  ! node-local rank comes from the launcher (OpenMPI or SLURM); without it
  ! the global rank is used, which is right on a single node.
  subroutine select_device()
    integer :: num_dev, dev, local_rank, length, status
    character(len=32) :: text
    num_dev = omp_get_num_devices()
    if (num_dev < 1) then
      print *, 'ERROR: rank', iproc, 'sees no OpenMP target device'
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end if
    local_rank = iproc
    call get_environment_variable('OMPI_COMM_WORLD_LOCAL_RANK', text, length, status)
    if (status /= 0) call get_environment_variable('SLURM_LOCALID', text, length, status)
    if (status == 0 .and. length > 0) read (text, *) local_rank
    dev = mod(local_rank, num_dev)
    call omp_set_default_device(dev)
    print '(A,I4,A,I4,A,I2,A,I2)', ' rank', iproc, ' local rank', local_rank, ' uses device', dev, ' of', num_dev
    !$omp target
    if (omp_is_initial_device()) print *, 'WARNING: target region ran on the host'
    !$omp end target
  end subroutine select_device
#endif

end program hst
