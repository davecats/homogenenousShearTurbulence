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
! Written and read collectively with MPI-IO: each rank repacks the rows
! ny0..nyN of its x slab into the CPL index order and writes them through
! a subarray view of the file (cpl_view); the file's four ghost rows are
! written by the ranks that own the rows they copy, in two more collective
! writes.  Reading takes the interior rows; the ghost rows are refilled by
! fill_ghosts.
module hst_io

  use, intrinsic :: iso_c_binding
  use mpi_f08
  use hst_params
  use hst_initial, only: generate_initial_field
  use hst_derivatives, only: s2_of, gamma_y_of, fill_ghosts

  implicit none
  private
  public :: restart_read, restart_write, field_write, make_output_dirs

  character, parameter :: LF = achar(10), TAB = achar(9)
  integer, parameter :: CPL_ORDER(3) = [1, 3, 2]     ! our component behind CPL's (u, v, w)

contains

  ! fields/ and p_fields/ next to the deck, as the CPL code expects; spectra/ for hst_stats.
  subroutine make_output_dirs()
    if (has_terminal) then
      call execute_command_line('mkdir -p fields p_fields spectra')
    end if
  end subroutine make_output_dirs

  ! Reads filename into V, or generates the initial field when it is absent.
  ! With time_from_restart the clock is taken from the file.  A restart may
  ! have a different resolution when restart_remap='mode_index': x and z are
  ! Fourier indices, so only modes present in both fields are copied; y is
  ! resampled by nearest neighbour in physical space.  The source CPL field
  ! is read only for the x slab needed by this rank.
  subroutine restart_read(filename)
    character(len=*), intent(in) :: filename
    integer :: io, ierr, unit, p
    character(len=4096) :: head
    integer(C_INT) :: r_nx, r_ny, r_nz, source_ny
    real(C_DOUBLE) :: r_alfa0, r_beta0, r_re, r_time, r_S
    integer(MPI_OFFSET_KIND) :: disp
    type(MPI_File) :: fh
    type(MPI_Datatype) :: view
    complex(C_DOUBLE_COMPLEX), allocatable :: buf(:, :, :, :)
    integer(C_INT) :: ix, iy, iz, src_iy
    integer :: sx0, sx1, sxB, source_zB
    real(C_DOUBLE) :: source_dy, target_y

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
    source_ny = r_nz - 1
    if (r_nx < 0 .or. r_ny < 0 .or. source_ny < 1) then
      if (has_terminal) print *, 'ERROR: invalid mesh dimensions in '//trim(filename)
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end if
    if (has_terminal .and. (r_nx /= nx .or. r_ny /= nz .or. source_ny /= ny)) then
      print *, '   interpolating restart mesh: file (nx, nz, ny) = ', r_nx, r_ny, source_ny
      print *, '                              deck (nx, nz, ny) = ', nx, nz, ny
    end if
    if (abs(r_alfa0 - alfa0) > 1.0d-4*max(abs(alfa0), 1.0d0) .or. &
        abs(r_beta0 - beta0) > 1.0d-4*max(abs(beta0), 1.0d0)) then
      if (trim(restart_remap) == 'mode_index') then
        if (has_terminal) print *, '   restart_remap=mode_index: remapping signed Fourier indices and projecting on target wavenumbers'
      else
      if (has_terminal) print *, 'ERROR: restart alpha0/beta0 do not match the deck:', r_alfa0, r_beta0, alfa0, beta0
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
      end if
    end if
    if (has_terminal .and. (abs(r_re - re) > 1.0d-9*re .or. r_S /= S)) &
      print *, '   note: Re or S differ from the file (', r_re, r_S, ')'
    if (time_from_restart) time = r_time

    ! The source file stores x=0..r_nx, z=-r_ny..r_ny and interior shear
    ! rows 0..source_ny-1 at file rows 2..r_nz.  Read the intersection of
    ! this rank's x slab and the source x range.  The other ranks still join
    ! the collective read with a zero-length operation.
    sx0 = max(int(nx0), 0)
    sx1 = min(int(nxN), int(r_nx))
    sxB = max(0, sx1 - sx0 + 1)
    source_zB = 2*int(r_ny) + 1
    call MPI_File_open(MPI_COMM_WORLD, trim(filename), MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr)
    if (sxB > 0) then
      allocate (buf(3, 0:source_ny-1, -r_ny:r_ny, sx0:sx1))
      call MPI_Type_create_subarray(4, [3, int(r_nz)+3, source_zB, int(r_nx)+1], &
                                    [3, int(source_ny), source_zB, sxB], [0, 2, 0, sx0], &
                                    MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, view, ierr)
      call MPI_Type_commit(view, ierr)
      call MPI_File_set_view(fh, disp, MPI_DOUBLE_COMPLEX, view, 'native', MPI_INFO_NULL, ierr)
      call MPI_File_read_all(fh, buf, size(buf), MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, ierr)
      call MPI_Type_free(view, ierr)
    else
      allocate (buf(3, 1, 1, 1))
      call MPI_File_set_view(fh, disp, MPI_DOUBLE_COMPLEX, MPI_DOUBLE_COMPLEX, 'native', MPI_INFO_NULL, ierr)
      call MPI_File_read_all(fh, buf, 0, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, ierr)
    end if
    call MPI_File_close(fh, ierr)

    ! The old CPL shear grid is uniform, irrespective of the target grid's
    ! optional stretching.  Use coordinates, not a ratio of array indices,
    ! and wrap the endpoint periodically before selecting the nearest row.
    source_dy = ly/real(source_ny, C_DOUBLE)
    do ix = nx0, nxN
      if (ix < sx0 .or. ix > sx1) cycle
      do iz = max(-nz, -r_ny), min(nz, r_ny)
        do iy = ny0, nyN
          target_y = modulo(y(iy), ly)
          src_iy = nint(target_y/source_dy)
          if (src_iy == source_ny) src_iy = 0
          V(iy, iz, ix, CPL_ORDER) = buf(:, src_iy, iz, ix)
        end do
      end do
    end do
    deallocate (buf)
    if (trim(restart_remap) == 'mode_index') then
      ! The host copy is now complete.  Ghost rows must be made with the
      ! target wavenumbers before the discrete divergence is evaluated.
      !$omp target update to(V)
      call fill_ghosts(1, 3)
      call fill_ghosts(2)
      call project_restart()
    end if
  end subroutine restart_read

  ! Enforce the target-grid discrete divergence constraint after a mode-index
  ! remap.  In the D0-weighted compact scheme the constraint is
  !   i*alfa*u + D1*v + i*beta*w = 0.
  ! For each mode we change only one component, avoiding a singular division
  ! at (alfa,beta)=(0,0).  This is intentionally a remap repair, not a general
  ! pressure projection: the source field is already close to divergence free.
  subroutine project_restart()
    integer(C_INT) :: ix, iz, iy, j
    complex(C_DOUBLE_COMPLEX) :: div, du, dv, dw
    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(V, der, ialfa, ibeta, nx0, nxN, nz, ny0, nyN) &
    !$omp private(ix, iz, iy, j, div, du, dv, dw)
    do ix = nx0, nxN
      do iz = -nz, nz
        do iy = ny0, nyN
          dv = 0.0d0
          du = 0.0d0
          dw = 0.0d0
          do j = -2, 2
            du = du + der(iy, 0, j)*V(iy + j, iz, ix, 1)
            dv = dv + der(iy, 1, j)*V(iy + j, iz, ix, 2)
            dw = dw + der(iy, 0, j)*V(iy + j, iz, ix, 3)
          end do
          if (ix /= 0) then
            div = dv + ibeta(iz)*dw + ialfa(ix)*du
            V(iy, iz, ix, 1) = V(iy, iz, ix, 1) - div/ialfa(ix)
          else if (iz /= 0) then
            div = dv + ibeta(iz)*dw
            V(iy, iz, ix, 3) = V(iy, iz, ix, 3) - div/ibeta(iz)
          else
            ! At the mean mode incompressibility only requires D1*v=0.
            V(iy, iz, ix, 2) = 0.0d0
          end if
        end do
      end do
    end do
    !$omp target update from(V)
    call fill_ghosts(1, 3)
    call fill_ghosts(2)
  end subroutine project_restart

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
    character(len=:), allocatable :: head
    complex(C_DOUBLE_COMPLEX), allocatable :: buf(:, :, :, :)
    integer(C_INT) :: ix, iy, iz

    ! header, built on the terminal rank (raw8 puts binary bytes in it: no trim)
    head = ''
    if (has_terminal) &
      head = 'nx='//str_i(nx)//' '//TAB//'ny='//str_i(nz)//' '//TAB//'nz='//str_i(ny + 1)// &
             ' '//TAB//'alpha0='//str_r(alfa0)//' '//TAB//'beta0='//str_r(beta0)// &
             ' '//TAB//'htcoeff='//str_r(merge(ystretch, -1.0d0, ystretch > 0.0d0))//' '//TAB//'Re='//str_r(re)//' '//TAB//'Pr=0.71'//LF// &
             'deltat='//str_r(deltat)//' '//TAB//'t_max='//str_r(t_max)//' '//TAB// &
             'dt_field='//str_r(dt_field)//' '//TAB//'dt_save='//str_r(dt_save)//LF// &
             't_field='//str_r(time)//LF// &
             'meanpx=0 '//TAB//'meanflowx=0 '//TAB//'meanpy=0 '//TAB//'meanflowy=0'//LF// &
             'time='//LF//raw8(time)//LF//'S='//LF//raw8(S)//LF//'S2='//LF//raw8(s2_of(time))//LF// &
             'gamma_x='//LF//raw8(S*time)//LF//'gamma_y='//LF//raw8(gamma_y_of(time))//LF//'Vfield='//LF
    ! this rank's rows in CPL order
    allocate (buf(3, ny0:nyN, -nz:nz, nx0:nxN))
    do ix = nx0, nxN
      do iz = -nz, nz
        do iy = ny0, nyN
          buf(:, iy, iz, ix) = V(iy, iz, ix, CPL_ORDER)
        end do
      end do
    end do
    call write_cpl(filename, head, buf, 3)
    deallocate (buf)
  end subroutine restart_write

  ! Writes one field with the layout of a component of V (host copy) as a
  ! headerless CPL array: the pressure files of prepare_pressure.cpl.
  subroutine field_write(filename, field)
    character(len=*), intent(in) :: filename
    complex(C_DOUBLE_COMPLEX), intent(in) :: field(ny0 - 2:, -nz:, nx0:)
    complex(C_DOUBLE_COMPLEX), allocatable :: buf(:, :, :, :)
    integer(C_INT) :: ix, iy, iz

    allocate (buf(1, ny0:nyN, -nz:nz, nx0:nxN))
    do ix = nx0, nxN
      do iz = -nz, nz
        do iy = ny0, nyN
          buf(1, iy, iz, ix) = field(iy, iz, ix)
        end do
      end do
    end do
    call write_cpl(filename, '', buf, 1)
    deallocate (buf)
  end subroutine field_write

  ! Creates filename, writes the header (built on the terminal rank; empty
  ! for none) and then buf, this rank's rows in CPL order with ncomp
  ! components: the interior rows, then the file's ghost rows ny, ny+1
  ! (copies of the rows 0, 1) and -2, -1 (copies of ny-2, ny-1) from the
  ! ranks that own those rows.
  subroutine write_cpl(filename, head, buf, ncomp)
    character(len=*), intent(in) :: filename, head
    integer, intent(in) :: ncomp
    complex(C_DOUBLE_COMPLEX), intent(in) :: buf(ncomp, ny0:nyN, -nz:nz, nx0:nxN)
    type(MPI_File) :: fh
    type(MPI_Status) :: status
    integer :: ierr, hlen
    integer(MPI_OFFSET_KIND) :: disp

    hlen = len(head)
    call MPI_Bcast(hlen, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    disp = hlen
    call MPI_File_open(MPI_COMM_WORLD, trim(filename), ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), MPI_INFO_NULL, fh, ierr)
    call MPI_File_set_size(fh, 0_MPI_OFFSET_KIND, ierr)
    if (has_terminal .and. hlen > 0) call MPI_File_write(fh, head, hlen, MPI_CHARACTER, status, ierr)
    call write_rows(fh, disp, buf, ncomp, ny0, nyN, ny0, .true.)
    call write_rows(fh, disp, buf, ncomp, 0, 1, ny, ny0 == 0)
    call write_rows(fh, disp, buf, ncomp, ny - 2, ny - 1, -2, nyN == ny - 1)
    call MPI_File_close(fh, ierr)
  end subroutine write_cpl

  ! Collective write of the rows b0..b1 of buf to the file rows f0.. (CPL
  ! row index -2..ny+1); a rank with nothing to write (active false) takes
  ! part with an empty write.
  subroutine write_rows(fh, disp, buf, ncomp, b0, b1, f0, active)
    type(MPI_File), intent(in) :: fh
    integer(MPI_OFFSET_KIND), intent(in) :: disp
    integer, intent(in) :: ncomp, b0, b1, f0
    complex(C_DOUBLE_COMPLEX), intent(in) :: buf(ncomp, ny0:nyN, -nz:nz, nx0:nxN)
    logical, intent(in) :: active
    type(MPI_Datatype) :: filetype, memtype
    type(MPI_Status) :: status
    integer :: ierr

    if (.not. active) then
      call MPI_File_set_view(fh, disp, MPI_DOUBLE_COMPLEX, MPI_DOUBLE_COMPLEX, 'native', MPI_INFO_NULL, ierr)
      call MPI_File_write_all(fh, buf, 0, MPI_DOUBLE_COMPLEX, status, ierr)
      return
    end if
    filetype = cpl_view(ncomp, f0, f0 + b1 - b0)
    call MPI_Type_create_subarray(4, [ncomp, nyB, 2*nz + 1, nxB], [ncomp, b1 - b0 + 1, 2*nz + 1, nxB], &
                                  [0, b0 - ny0, 0, 0], MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, memtype, ierr)
    call MPI_Type_commit(memtype, ierr)
    call MPI_File_set_view(fh, disp, MPI_DOUBLE_COMPLEX, filetype, 'native', MPI_INFO_NULL, ierr)
    call MPI_File_write_all(fh, buf, 1, memtype, status, ierr)
    call MPI_Type_free(memtype, ierr)
    call MPI_Type_free(filetype, ierr)
  end subroutine write_rows

  ! MPI-IO view of the file rows r0..r1 (CPL row index -2..ny+1) of this
  ! rank's x slab: the file array is (ncomp, ny+4, 2nz+1, nx+1) in Fortran
  ! order (ncomp = 3 for the velocity, 1 for the pressure).  Committed; the
  ! caller frees it.
  function cpl_view(ncomp, r0, r1) result(view)
    integer, intent(in) :: ncomp, r0, r1
    type(MPI_Datatype) :: view
    integer :: ierr
    call MPI_Type_create_subarray(4, [ncomp, ny + 4, 2*nz + 1, nx + 1], [ncomp, r1 - r0 + 1, 2*nz + 1, nxB], &
                                  [0, r0 + 2, 0, nx0], MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, view, ierr)
    call MPI_Type_commit(view, ierr)
  end function cpl_view

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
