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
  public :: T_TRANSFORM, T_PREPARE, T_PRODUCTS, T_SHIFT, T_SOLVE, T_RECOVER, T_OTHER

  integer, parameter :: T_TRANSFORM = 1, T_PREPARE = 2, T_PRODUCTS = 3, T_SHIFT = 4, T_SOLVE = 5, &
                        T_RECOVER = 6, T_OTHER = 7, NPHASE = 7
  character(len=24), parameter :: names(NPHASE) = [character(len=24) :: &
    'transform to physical', 'buildrhs_prepare', 'products + buildrhs', 'shear_shift', &
    'implicit solves', 'ghosts, dv/dy, u and w', 'statistics, I/O, rest']
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
