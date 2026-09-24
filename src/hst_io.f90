! Restart files and field snapshots.
!
! One binary stream file: a header of 3 integers and 6 doubles
!   nx ny nz  alfa0 beta0 ly re S time
! followed by the complex array (ny, 2*nz+1, nx+1, 3) of (u, v, w) modes,
! rows 0..ny-1 only (the ghost rows are images and are rebuilt after a
! read).  Written and read collectively with MPI-IO.  A pressure snapshot
! (field_write) has the same header and one component, (ny, 2*nz+1, nx+1).
!
! From channel/src/io/restart_io.f90.
module hst_io

  use, intrinsic :: iso_c_binding
  use mpi_f08
  use hst_params
  use hst_mpi, only: file_view_type, memory_type, file_view_type1, memory_type1
  use hst_initial, only: generate_initial_field

  implicit none
  private
  public :: restart_read, restart_write, field_write

  integer(MPI_OFFSET_KIND), parameter :: header_bytes = 3*C_INT + 6*C_DOUBLE

contains

  ! Reads filename into V, or generates the initial field when it is absent.
  ! With time_from_restart the clock is taken from the file.
  subroutine restart_read(filename)
    character(len=*), intent(in) :: filename
    integer :: io, ierr
    integer(C_INT) :: r_nx, r_ny, r_nz
    real(C_DOUBLE) :: r_alfa0, r_beta0, r_ly, r_re, r_S, r_time
    type(MPI_File) :: fh

    open (unit=120, file=trim(filename), access='stream', status='old', action='read', iostat=io)
    if (io /= 0) then
      if (has_terminal) print *, '   restart file '//trim(filename)//' not found'
      call generate_initial_field()
      return
    end if
    if (has_terminal) print *, '   reading '//trim(filename)
    read (120, pos=1) r_nx, r_ny, r_nz, r_alfa0, r_beta0, r_ly, r_re, r_S, r_time
    close (120)
    if (r_nx /= nx .or. r_ny /= ny .or. r_nz /= nz .or. r_alfa0 /= alfa0 .or. r_beta0 /= beta0 .or. r_ly /= ly) then
      if (has_terminal) then
        print *, 'ERROR: mesh in '//trim(filename)//' does not match the deck'
        print *, '   file: ', r_nx, r_ny, r_nz, r_alfa0, r_beta0, r_ly
        print *, '   deck: ', nx, ny, nz, alfa0, beta0, ly
      end if
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end if
    if (has_terminal .and. (r_re /= re .or. r_S /= S)) &
      print *, '   note: re or S differ from the restart file (', r_re, r_S, ')'
    if (time_from_restart) time = r_time

    call MPI_File_open(MPI_COMM_WORLD, trim(filename), MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr)
    call MPI_File_set_view(fh, header_bytes, MPI_DOUBLE_COMPLEX, file_view_type, 'native', MPI_INFO_NULL, ierr)
    call MPI_File_read_all(fh, V, 1, memory_type, MPI_STATUS_IGNORE, ierr)
    call MPI_File_close(fh, ierr)
  end subroutine restart_read

  ! Writes V (host copy) with the header.  The caller updates V from the
  ! device first.
  subroutine restart_write(filename)
    character(len=*), intent(in) :: filename
    type(MPI_File) :: fh
    type(MPI_Status) :: status
    integer :: ierr

    call MPI_File_open(MPI_COMM_WORLD, trim(filename), ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), MPI_INFO_NULL, fh, ierr)
    if (has_terminal) then
      call MPI_File_write(fh, [nx, ny, nz], 3, MPI_INTEGER, status, ierr)
      call MPI_File_write(fh, [alfa0, beta0, ly, re, S, time], 6, MPI_DOUBLE_PRECISION, status, ierr)
    end if
    call MPI_File_set_view(fh, header_bytes, MPI_DOUBLE_COMPLEX, file_view_type, 'native', MPI_INFO_NULL, ierr)
    call MPI_File_write_all(fh, V, 1, memory_type, status, ierr)
    call MPI_File_close(fh, ierr)
  end subroutine restart_write

  ! Writes one field with the layout of a component of V (host copy), with
  ! the same header as a restart file.
  subroutine field_write(filename, field)
    character(len=*), intent(in) :: filename
    complex(C_DOUBLE_COMPLEX), intent(in) :: field(ny0 - 2:, -nz:, nx0:)
    type(MPI_File) :: fh
    type(MPI_Status) :: status
    integer :: ierr

    call MPI_File_open(MPI_COMM_WORLD, trim(filename), ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), MPI_INFO_NULL, fh, ierr)
    if (has_terminal) then
      call MPI_File_write(fh, [nx, ny, nz], 3, MPI_INTEGER, status, ierr)
      call MPI_File_write(fh, [alfa0, beta0, ly, re, S, time], 6, MPI_DOUBLE_PRECISION, status, ierr)
    end if
    call MPI_File_set_view(fh, header_bytes, MPI_DOUBLE_COMPLEX, file_view_type1, 'native', MPI_INFO_NULL, ierr)
    call MPI_File_write_all(fh, field, 1, memory_type1, status, ierr)
    call MPI_File_close(fh, ierr)
  end subroutine field_write

end module hst_io
