! Wall-clock time per phase of the substep, accumulated over the run and
! printed at the end when timing = .true. in &time_control.  Each toc
! waits for the device first, so the phases are what the GPU spent, not
! what the host spent launching.  With timing off the calls return at once.
!
!   call tic()          start (or restart) the clock
!   call toc(T_X)       add the time since the last tic/toc to phase T_X
!   call timer_report(nsteps)
module hst_timer

  use, intrinsic :: iso_c_binding
  use mpi_f08, only: MPI_Wtime
  use hst_params, only: timing, has_terminal
  use hst_fft, only: device_sync

  implicit none
  private
  public :: tic, toc, timer_report
  public :: T_TRANSFORM, T_PACK, T_ALLTOALL, T_PREPARE, T_PRODUCTS, T_SHIFT, T_SOLVE, T_RECOVER, T_OTHER

  ! The two transpose phases are charged from inside hst_mpi, so the
  ! transform and product phases are the FFTs and kernels around them.  On
  ! one rank "pack, unpack" is the local repack and there is no alltoall.
  integer, parameter :: T_TRANSFORM = 1, T_PACK = 2, T_ALLTOALL = 3, T_PREPARE = 4, T_PRODUCTS = 5, &
                        T_SHIFT = 6, T_SOLVE = 7, T_RECOVER = 8, T_OTHER = 9, NPHASE = 9
  character(len=24), parameter :: names(NPHASE) = [character(len=24) :: &
    'to physical: FFTs, CFL', 'transpose pack, unpack', 'transpose alltoall', 'buildrhs_prepare', &
    'products, FFTs, buildrhs', 'shear_shift', 'implicit solves', 'ghosts, dv/dy, u and w', &
    'statistics, I/O, rest']
  real(C_DOUBLE), save :: acc(NPHASE) = 0.0d0, t_last = 0.0d0

contains

  subroutine tic()
    if (.not. timing) return
    call device_sync()
    t_last = MPI_Wtime()
  end subroutine tic

  subroutine toc(phase)
    integer, intent(in) :: phase
    real(C_DOUBLE) :: t
    if (.not. timing) return
    call device_sync()
    t = MPI_Wtime()
    acc(phase) = acc(phase) + (t - t_last)
    t_last = t
  end subroutine toc

  ! Per-phase totals of this rank (rank 0 prints; ranks differ by the wait
  ! inside the alltoall only), per step and as a share of the sum.
  subroutine timer_report(nsteps)
    integer(C_SIZE_T), intent(in) :: nsteps
    integer :: i
    real(C_DOUBLE) :: total
    if (.not. (timing .and. has_terminal)) return
    total = sum(acc)
    write (*, '(A)') ' '
    write (*, '(A,I0,A)') '   time per phase over ', nsteps, ' steps (rank 0):'
    do i = 1, NPHASE
      write (*, '(A,A,F10.3,A,F9.5,A,F6.1,A)') '     ', names(i), acc(i), ' s', acc(i)/max(nsteps, 1_C_SIZE_T), &
        ' s/step', 100.0d0*acc(i)/max(total, tiny(total)), ' %'
    end do
    write (*, '(A,A,F10.3,A,F9.5,A)') '     ', 'total                   ', total, ' s', total/max(nsteps, 1_C_SIZE_T), ' s/step'
  end subroutine timer_report

end module hst_timer
