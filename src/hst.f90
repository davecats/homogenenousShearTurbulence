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
! Set-up, then per step: timestep() (hst_equations), statistics, snapshots,
! restart file, new time step from the CFL number.
!
program hst

  use, intrinsic :: iso_c_binding
  use mpi_f08
  use hst_params
  use hst_input
  use hst_mpi
  use hst_fft
  use hst_setup
  use hst_derivatives
  use hst_linsolve
  use hst_transforms
  use hst_equations
  use hst_io
  use hst_stats
  use hst_pressure
  use hst_stokes
#ifdef HAVE_CUDA
  use omp_lib
#endif

  implicit none

  character(len=256) :: deck
  character(len=40) :: fname
  integer :: ierr, m
  real(C_DOUBLE) :: cfl_global, t0, t1, elapsed

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, iproc, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
  has_terminal = (iproc == 0)
#ifdef HAVE_CUDA
  call select_device()
#endif

  !------------------------------------------------------------ set-up ----
  deck = 'hst.in'
  if (command_argument_count() >= 1) call get_command_argument(1, deck)
  call read_input(trim(deck))
  if (has_terminal) call print_input()
  call setup_decomposition()
  call allocate_fields()
  call init_fft()
  call init_linsolve()
  call setup_derivatives()
  call make_output_dirs()
  call restart_read('Dati.cart.out')
  !$omp target update to(V)
  call stokes_setup()
  call stokes_apply()
  do m = 1, 3
    call fill_ghosts(m)
  end do
  call open_runtimedata()

  ! First time step from the CFL of the initial field.
  if (deltat == 0.0d0) deltat = 1.0d0
  call transform_to_physical()
  call compute_cfl()
  call new_timestep()
  ifield = floor(time/dt_field)          ! CPL numbering: fields/field<ifield+1>.fld is the next one
  if (has_terminal) write (*, '(A)') '        time       deltat       cfl         energy           diss' // &
    '           uw/2           vw/2      (CPL Runtimedata columns 10-13)'
  call outstats()

  !--------------------------------------------------------- time loop ----
  elapsed = 0.0d0
  do while (time < t_max - 0.5d0*deltat .and. istep < nstep)
    t0 = MPI_Wtime()
    istep = istep + 1
    call timestep()

    if (crossed(dt_stat)) call outstats()
    if (crossed(dt_field)) then
      ifield = ifield + 1
      write (fname, '(A,I0,A)') 'fields/field', ifield, '.fld'
      if (has_terminal) print '(A,F12.5)', '   writing '//trim(fname)//' at time', time
      !$omp target update from(V)
      call restart_write(trim(fname))
      write (fname, '(A,I0,A)') 'p_fields/pField', ifield, '.fld'
      call write_pressure(trim(fname))
    end if
    if (crossed(dt_save)) then
      if (has_terminal) print '(A,F12.5)', '   writing Dati.cart.out at time', time
      !$omp target update from(V)
      call restart_write('Dati.cart.out')
    end if
    call new_timestep()

    t1 = MPI_Wtime()
    elapsed = elapsed + (t1 - t0)
    if (has_terminal .and. (mod(istep, 50_C_SIZE_T) == 0 .or. istep <= 5)) &
      write (*, '(A,I0,A,F9.5,A,F12.2,A)') '   step ', istep, ': ', t1 - t0, ' s/step, ', elapsed, ' s elapsed'
  end do

  !------------------------------------------------------------ finish ----
  if (has_terminal) print '(A,F12.5,A,I0,A)', '   end of run at time', time, ' after ', istep, ' steps; writing Dati.cart.out'
  !$omp target update from(V)
  call restart_write('Dati.cart.out')
  call close_runtimedata()
  call free_linsolve()
  call free_fft()
  call free_fields()
  call free_mpi()
  call MPI_Finalize(ierr)

contains

  ! .true. when the interval boundary of dt is crossed by this step.
  logical function crossed(dt)
    real(C_DOUBLE), intent(in) :: dt
    crossed = .false.
    if (dt > 0.0d0) crossed = floor((time + 0.5d0*deltat)/dt) > floor((time - 0.5d0*deltat)/dt)
  end function crossed

  ! Reduce the CFL number over ranks and choose the next time step: from
  ! cflmax when the deck gives deltat = 0, else the fixed value (capped by
  ! cflmax if that is set).
  subroutine new_timestep()
    call MPI_Allreduce(cfl, cfl_global, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    cfl = cfl_global
    if (cflmax > 0.0d0 .and. cfl > 0.0d0) deltat = cflmax/cfl
    if (dt_fixed > 0.0d0) deltat = min(deltat, dt_fixed)
    if (dt_fixed > 0.0d0 .and. cflmax <= 0.0d0) deltat = dt_fixed
  end subroutine new_timestep

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
