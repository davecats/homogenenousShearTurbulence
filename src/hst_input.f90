! The input deck: one Fortran namelist file, four groups.  Every rank reads
! it.  Defaults are set here, so a deck only has to name what it changes.
module hst_input

  use, intrinsic :: iso_c_binding
  use hst_params

  implicit none
  private
  public :: read_input, print_input, fft_fit

contains

  subroutine read_input(filename)
    character(len=*), intent(in) :: filename
    integer :: unit, ios
    real(C_DOUBLE), parameter :: PI = 3.141592653589793d0

    namelist /mesh/ nx, ny, nz, alfa0, beta0, ly
    namelist /physics/ re, S, linear, exact_shift
    namelist /time_control/ deltat, cflmax, t_max, nstep, dt_stat, dt_field, dt_save, time, time_from_restart
    namelist /init/ amplitude, seed, kpeak

    ! defaults: the Sekimoto, Dong & Jimenez (2016) box Lx:Ly:Lz = 3:2:1
    nx = 63; ny = 128; nz = 63
    ly = 2.0d0; alfa0 = 2.0d0*PI/3.0d0; beta0 = 2.0d0*PI
    re = 1000.0d0; S = 1.0d0; linear = .false.; exact_shift = .false.
    deltat = 0.0d0; cflmax = 1.0d0; t_max = 100.0d0; nstep = 1000000
    dt_stat = 0.01d0; dt_field = 10.0d0; dt_save = 10.0d0
    time = 0.0d0; time_from_restart = .false.
    amplitude = 1.0d-3; seed = 1; kpeak = 4.0d0

    open (newunit=unit, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      print *, 'ERROR: cannot open input file '//trim(filename)
      error stop 1
    end if
    read (unit, nml=mesh); rewind (unit)
    read (unit, nml=physics); rewind (unit)
    read (unit, nml=time_control); rewind (unit)
    read (unit, nml=init)
    close (unit)

    ! derived quantities
    ni = 1.0d0/re
    dt_fixed = deltat
    lx = 2.0d0*PI/alfa0
    lz = 2.0d0*PI/beta0
    nxd = 3*(nx + 1)/2; do while (.not. fft_fit(nxd)); nxd = nxd + 1; end do
    nzd = 3*nz;         do while (.not. fft_fit(nzd)); nzd = nzd + 1; end do
    dx = PI/(alfa0*nxd)
    dz = 2.0d0*PI/(beta0*nzd)
    dy = ly/ny
    factor = 1.0d0/(2.0d0*nxd*nzd)
  end subroutine read_input

  ! .true. when n is a power of two times at most one factor three, which
  ! is what keeps the FFTs fast.
  logical function fft_fit(n)
    integer(C_INT), intent(in) :: n
    integer(C_INT) :: j
    j = n
    do while (mod(j, 2) == 0)
      j = j/2
    end do
    fft_fit = (j == 1 .or. j == 3)
  end function fft_fit

  subroutine print_input()
    write (*, '(A)') ' '
    write (*, '(A)') '!====================================================!'
    write (*, '(A)') '!     DNS of homogeneous shear turbulence  (hst)     !'
    write (*, '(A)') '!====================================================!'
    write (*, '(A)') ' '
    write (*, '(A,I6,A,I6,A,I6)') '   nx    =', nx, '   ny    =', ny, '   nz    =', nz
    write (*, '(A,I6,A,I6)') '   nxd   =', nxd, '   nzd   =', nzd
    write (*, '(A,F10.6,A,F10.6,A,F10.6)') '   lx    =', lx, '   ly    =', ly, '   lz    =', lz
    write (*, '(A,F10.6,A,F10.6)') '   alfa0 =', alfa0, '   beta0 =', beta0
    write (*, '(A,F10.2,A,F10.6,A,L1,A,L1)') '   re    =', re, '   S     =', S, '   linear =', linear, &
      '   exact_shift =', exact_shift
    write (*, '(A,F10.6,A,F10.6,A,F10.3)') '   deltat=', deltat, '   cflmax=', cflmax, '   t_max =', t_max
    write (*, '(A,F10.4,A,F10.4,A,F10.4)') '   dt_stat=', dt_stat, '  dt_field=', dt_field, '  dt_save=', dt_save
    write (*, '(A,I10,A,L1)') '   nstep =', nstep, '   time_from_restart = ', time_from_restart
    write (*, '(A,ES10.3,A,I6,A,F8.3)') '   amplitude=', amplitude, '   seed =', seed, '   kpeak =', kpeak
    write (*, '(A)') ' '
  end subroutine print_input

end module hst_input
