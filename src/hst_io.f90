! Field files in the layout of the CPL code hst-main, so that its
! post-processing chain (postprocess/, pressure_reconstruction/) reads them
! unchanged and its fields can be read here.
!
! Velocity (restart file Dati.cart.out and snapshots fields/field<n>.fld):
! a text header exactly as io.cpl writes it,
!   nx=.. \tny=.. \tnz=.. \talpha0=.. \tbeta0=.. \thtcoeff=-1 \tRe=.. \tPr=0.71
!   deltat=.. \tt_max=.. \tdt_field=.. \tdt_save=..
!   t_field=..
!   meanpx=0 \tmeanflowx=0 \tmeanpy=0 \tmeanflowy=0
!   time=      <8 raw bytes>
!   S=         <8 raw bytes>
!   S2=        <8 raw bytes>
!   gamma_x=   <8 raw bytes>
!   gamma_y=   <8 raw bytes>
!   Vfield=
! then the array  ARRAY(0..nx, -ny_cpl..ny_cpl, -1..nz_cpl+1) OF (u, v, w)
! complex, stored with the last index fastest (C order) and the three
! components innermost.  CPL names: ny_cpl = our nz (spanwise modes),
! nz_cpl = our ny + 1 (their nz-1 unique points over 2 = our ny points over
! ly), their (v, w) = our (w, v), and their row iz_cpl is our iy = iz_cpl - 1.
! The four ghost rows are written as plain periodic copies (the CPL code
! re-applies its periodic condition at start-up; ours refills them).
!
! Pressure (p_fields/pField<n>.fld): the same array without header and with
! one component, as prepare_pressure.cpl writes it.
!
! Written and read collectively with MPI-IO: each rank repacks its x slab
! into the CPL index order and writes it through a subarray view.
module hst_io

  use, intrinsic :: iso_c_binding
  use mpi_f08
  use hst_params
  use hst_mpi, only: cpl_view_type, cpl_pview_type
  use hst_initial, only: generate_initial_field
  use hst_derivatives, only: s2_of, gamma_y_of

  implicit none
  private
  public :: restart_read, restart_write, field_write, make_output_dirs

  character, parameter :: LF = achar(10), TAB = achar(9)

contains

  ! fields/ and p_fields/ next to the deck, as the CPL code expects.
  subroutine make_output_dirs()
    if (has_terminal) then
      call execute_command_line('mkdir -p fields p_fields')
    end if
  end subroutine make_output_dirs

  ! Reads filename into V, or generates the initial field when it is absent.
  ! With time_from_restart the clock is taken from the file.
  subroutine restart_read(filename)
    character(len=*), intent(in) :: filename
    integer :: io, ierr, unit, p
    character(len=4096) :: head
    integer(C_INT) :: r_nx, r_ny, r_nz
    real(C_DOUBLE) :: r_alfa0, r_beta0, r_re, r_time, r_S
    integer(MPI_OFFSET_KIND) :: disp
    type(MPI_File) :: fh
    complex(C_DOUBLE_COMPLEX), allocatable :: buf(:, :, :, :)
    integer(C_INT) :: ix, iy, iz

    open (newunit=unit, file=trim(filename), access='stream', status='old', action='read', iostat=io)
    if (io /= 0) then
      if (has_terminal) print *, '   restart file '//trim(filename)//' not found'
      call generate_initial_field()
      return
    end if
    if (has_terminal) print *, '   reading '//trim(filename)
    head = ''
    read (unit, pos=1, iostat=io) head
    p = index(head, 'Vfield='//LF)
    if (p == 0) then
      if (has_terminal) print *, 'ERROR: '//trim(filename)//' has no Vfield= line (not a CPL field file)'
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end if
    disp = p - 1 + 8
    r_nx = int(header_value(head, 'nx='))
    r_ny = int(header_value(head, 'ny='))
    r_nz = int(header_value(head, 'nz='))
    r_alfa0 = header_value(head, 'alpha0=')
    r_beta0 = header_value(head, 'beta0=')
    r_re = header_value(head, 'Re=')
    p = index(head, 'time='//LF); read (unit, pos=p + 6) r_time
    p = index(head, LF//'S='//LF); read (unit, pos=p + 4) r_S
    close (unit)
    if (r_nx /= nx .or. r_ny /= nz .or. r_nz /= ny + 1) then
      if (has_terminal) then
        print *, 'ERROR: mesh in '//trim(filename)//' does not match the deck'
        print *, '   file (CPL names nx ny nz): ', r_nx, r_ny, r_nz
        print *, '   deck (nx, nz, ny+1):       ', nx, nz, ny + 1
      end if
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end if
    if (has_terminal .and. (abs(r_alfa0 - alfa0) > 1.0d-4*alfa0 .or. abs(r_beta0 - beta0) > 1.0d-4*beta0)) &
      print *, '   note: alpha0/beta0 in the file differ from the deck:', r_alfa0, r_beta0
    if (has_terminal .and. (abs(r_re - re) > 1.0d-9*re .or. r_S /= S)) &
      print *, '   note: Re or S differ from the file (', r_re, r_S, ')'
    if (time_from_restart) time = r_time

    allocate (buf(3, ny0 - 2:nyN + 2, -nz:nz, nx0:nxN))
    call MPI_File_open(MPI_COMM_WORLD, trim(filename), MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr)
    call MPI_File_set_view(fh, disp, MPI_DOUBLE_COMPLEX, cpl_view_type, 'native', MPI_INFO_NULL, ierr)
    call MPI_File_read_all(fh, buf, size(buf), MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, ierr)
    call MPI_File_close(fh, ierr)
    do ix = nx0, nxN
      do iz = -nz, nz
        do iy = 0, ny - 1
          V(iy, iz, ix, 1) = buf(1, iy, iz, ix)
          V(iy, iz, ix, 3) = buf(2, iy, iz, ix)
          V(iy, iz, ix, 2) = buf(3, iy, iz, ix)
        end do
      end do
    end do
    deallocate (buf)
  end subroutine restart_read

  ! Number after `key` in the text header.
  real(C_DOUBLE) function header_value(head, key)
    character(len=*), intent(in) :: head, key
    integer :: p, q
    p = index(head, key)
    if (p == 0) then
      header_value = -huge(1.0d0)
      return
    end if
    p = p + len(key)
    q = p
    do while (q <= len(head) .and. head(q:q) /= ' ' .and. head(q:q) /= TAB .and. head(q:q) /= LF)
      q = q + 1
    end do
    read (head(p:q - 1), *) header_value
  end function header_value

  ! Writes V (host copy) with the CPL header.  The caller updates V from the
  ! device first.
  subroutine restart_write(filename)
    character(len=*), intent(in) :: filename
    type(MPI_File) :: fh
    type(MPI_Status) :: status
    integer :: ierr, hlen
    integer(MPI_OFFSET_KIND) :: disp
    character(len=1024) :: head
    complex(C_DOUBLE_COMPLEX), allocatable :: buf(:, :, :, :)
    integer(C_INT) :: ix, iy, iz, jy

    ! header, built on the terminal rank and written by it
    if (has_terminal) then
      head = 'nx='//str_i(nx)//' '//TAB//'ny='//str_i(nz)//' '//TAB//'nz='//str_i(ny + 1)// &
             ' '//TAB//'alpha0='//str_r(alfa0)//' '//TAB//'beta0='//str_r(beta0)// &
             ' '//TAB//'htcoeff='//str_r(merge(ystretch, -1.0d0, ystretch > 0.0d0))//' '//TAB//'Re='//str_r(re)//' '//TAB//'Pr=0.71'//LF// &
             'deltat='//str_r(deltat)//' '//TAB//'t_max='//str_r(t_max)//' '//TAB// &
             'dt_field='//str_r(dt_field)//' '//TAB//'dt_save='//str_r(dt_save)//LF// &
             't_field='//str_r(time)//LF// &
             'meanpx=0 '//TAB//'meanflowx=0 '//TAB//'meanpy=0 '//TAB//'meanflowy=0'//LF// &
             'time='//LF//raw8(time)//LF//'S='//LF//raw8(S)//LF//'S2='//LF//raw8(s2_of(time))//LF// &
             'gamma_x='//LF//raw8(S*time)//LF//'gamma_y='//LF//raw8(gamma_y_of(time))//LF//'Vfield='//LF
      hlen = index(head, 'Vfield='//LF) + 7
    end if
    call MPI_Bcast(hlen, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    disp = hlen

    ! this rank's slab in CPL order, with periodic ghost rows
    allocate (buf(3, ny0 - 2:nyN + 2, -nz:nz, nx0:nxN))
    do ix = nx0, nxN
      do iz = -nz, nz
        do iy = -2, ny + 1
          jy = modulo(iy, ny)
          buf(1, iy, iz, ix) = V(jy, iz, ix, 1)
          buf(2, iy, iz, ix) = V(jy, iz, ix, 3)
          buf(3, iy, iz, ix) = V(jy, iz, ix, 2)
        end do
      end do
    end do

    call MPI_File_open(MPI_COMM_WORLD, trim(filename), ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), MPI_INFO_NULL, fh, ierr)
    call MPI_File_set_size(fh, 0_MPI_OFFSET_KIND, ierr)
    if (has_terminal) call MPI_File_write(fh, head(1:hlen), hlen, MPI_CHARACTER, status, ierr)
    call MPI_File_set_view(fh, disp, MPI_DOUBLE_COMPLEX, cpl_view_type, 'native', MPI_INFO_NULL, ierr)
    call MPI_File_write_all(fh, buf, size(buf), MPI_DOUBLE_COMPLEX, status, ierr)
    call MPI_File_close(fh, ierr)
    deallocate (buf)
  end subroutine restart_write

  ! Writes one field with the layout of a component of V (host copy) as a
  ! headerless CPL array: the pressure files of prepare_pressure.cpl.
  subroutine field_write(filename, field)
    character(len=*), intent(in) :: filename
    complex(C_DOUBLE_COMPLEX), intent(in) :: field(ny0 - 2:, -nz:, nx0:)
    type(MPI_File) :: fh
    type(MPI_Status) :: status
    integer :: ierr
    complex(C_DOUBLE_COMPLEX), allocatable :: buf(:, :, :)
    integer(C_INT) :: ix, iy, iz

    allocate (buf(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN))
    do ix = nx0, nxN
      do iz = -nz, nz
        do iy = -2, ny + 1
          buf(iy, iz, ix) = field(modulo(iy, ny), iz, ix)
        end do
      end do
    end do
    call MPI_File_open(MPI_COMM_WORLD, trim(filename), ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), MPI_INFO_NULL, fh, ierr)
    call MPI_File_set_size(fh, 0_MPI_OFFSET_KIND, ierr)
    call MPI_File_set_view(fh, 0_MPI_OFFSET_KIND, MPI_DOUBLE_COMPLEX, cpl_pview_type, 'native', MPI_INFO_NULL, ierr)
    call MPI_File_write_all(fh, buf, size(buf), MPI_DOUBLE_COMPLEX, status, ierr)
    call MPI_File_close(fh, ierr)
    deallocate (buf)
  end subroutine field_write

  function str_i(i) result(s)
    integer(C_INT), intent(in) :: i
    character(len=:), allocatable :: s
    character(len=32) :: t
    write (t, '(I0)') i
    s = trim(t)
  end function str_i

  function str_r(x) result(s)
    real(C_DOUBLE), intent(in) :: x
    character(len=:), allocatable :: s
    character(len=32) :: t
    write (t, '(G0.17)') x
    s = trim(adjustl(t))
  end function str_r

  ! The 8 bytes of a double, as CPL's WRITE BINARY puts them in the header.
  function raw8(x) result(s)
    real(C_DOUBLE), intent(in) :: x
    character(len=8) :: s
    s = transfer(x, s)
  end function raw8

end module hst_io
