! What every test program does before and after its own work.
!
!   call test_start()          MPI, one device per rank, the deck (first
!                              argument, default hst.in); parameters may be
!                              changed after this call
!   call test_setup()          decomposition, fields, transforms, line
!                              solver, stencils
!   call test_finish(passed)   frees everything, prints PASSED or FAILED,
!                              exits with status 1 on failure
module test_common

  use, intrinsic :: iso_c_binding
  use mpi_f08
  use hst_params
  use hst_input
  use hst_mpi
  use hst_fft
  use hst_setup
  use hst_derivatives
  use hst_linsolve

  implicit none
  private
  public :: test_start, test_setup, test_finish

contains

  subroutine test_start()
    character(len=256) :: deck
    integer :: ierr
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, iproc, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
    has_terminal = (iproc == 0)
    call select_device()
    deck = 'hst.in'
    if (command_argument_count() >= 1) call get_command_argument(1, deck)
    call read_input(trim(deck))
  end subroutine test_start

  subroutine test_setup()
    call setup_decomposition()
    call allocate_fields()
    call init_fft()
    call init_linsolve()
    call setup_derivatives()
  end subroutine test_setup

  subroutine test_finish(passed)
    logical, intent(in) :: passed
    integer :: ierr
    call free_linsolve(); call free_fft(); call free_fields(); call free_mpi()
    call MPI_Finalize(ierr)
    if (.not. passed) then
      print *, 'FAILED'
      error stop 1
    end if
    if (has_terminal) print *, 'PASSED'
  end subroutine test_finish

end module test_common
