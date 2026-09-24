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
  use hst_mpi
  use hst_fft
  use hst_setup
  use hst_transforms
  use hst_io
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
  call allocate_fields()
  call init_fft()
  call restart_read('Dati.cart.out')
  !$omp target update to(V)

  ! WP2 and later: derivatives, ghost rows, time loop, statistics.

  call free_fft()
  call free_fields()
  call free_mpi()
  call MPI_Finalize(ierr)

contains

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
