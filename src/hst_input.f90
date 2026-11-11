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

    namelist /mesh/ nx, ny, nz, alfa0, beta0, ly, ystretch, line_chunk, transport, npy
    namelist /physics/ re, S, linear, exact_shift, s2_mean, s2_amplitude, s2_period, s2_start, &
      sl_amplitude, sl_period, sl_delta, sl_start, sl_bodyforce, sl_ramp
    namelist /time_control/ deltat, cflmax, t_max, nstep, wall_max, dt_stat, dt_field, dt_save, time, time_from_restart, restart_from_file, restart_remap, restart_file, timing
    namelist /init/ initial_condition, amplitude, seed, kpeak, tg_amplitude, tg_mode
    namelist /damping/ damp_kpar, damp_kperp, damp_rate, damp_hard, damp_edge, damp_norm, damp_mean, damp_follow, damp_angle

    ! defaults: the Sekimoto, Dong & Jimenez (2016) box Lx:Ly:Lz = 3:2:1
    nx = 63; ny = 128; nz = 63
    ly = 2.0d0; alfa0 = 2.0d0*PI/3.0d0; beta0 = 2.0d0*PI; ystretch = 0.0d0; line_chunk = 0; transport = 'auto'; npy = 0
    re = 1000.0d0; S = 1.0d0; linear = .false.; exact_shift = .false.
    s2_mean = 0.0d0; s2_amplitude = 0.0d0; s2_period = 0.0d0; s2_start = 0.0d0
    sl_amplitude = 0.0d0; sl_period = 1.0d0; sl_delta = 0.02d0; sl_start = 0.0d0; sl_bodyforce = .true.; sl_ramp = .false.
    deltat = 0.0d0; cflmax = 1.0d0; t_max = 100.0d0; nstep = 1000000; wall_max = 0.0d0
    dt_stat = 0.01d0; dt_field = 10.0d0; dt_save = 10.0d0
    time = 0.0d0; time_from_restart = .false.; timing = .false.
    restart_from_file = .true.          ! preserve default behaviour: always try to read restart_file
    restart_remap = 'strict'
    restart_file = 'Dati.cart.out'
    initial_condition = 'random'
    amplitude = 1.0d-3; seed = 1; kpeak = 4.0d0
    damp_kpar = 0.0d0; damp_kperp = 0.0d0; damp_rate = 0.0d0; damp_hard = .false.; damp_edge = 0.0d0
    damp_norm = 2.0d0; damp_mean = .false.; damp_follow = .true.; damp_angle = 0.0d0
    tg_amplitude = 1.0d0; tg_mode = 1

    open (newunit=unit, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      print *, 'ERROR: cannot open input file '//trim(filename)
      error stop 1
    end if
    read (unit, nml=mesh); rewind (unit)
    read (unit, nml=physics); rewind (unit)
    read (unit, nml=time_control); rewind (unit)
    read (unit, nml=init); rewind (unit)
    read (unit, nml=damping, iostat=ios)        ! optional group: absent (ios < 0) keeps the defaults
    if (ios > 0) then
      print *, 'ERROR: cannot read the &damping group of '//trim(filename)
      error stop 1
    end if
    close (unit)
    if (transport /= 'auto' .and. transport /= 'mpi' .and. transport /= 'nccl') then
      print *, 'ERROR: transport must be auto, mpi or nccl'
      error stop 1
    end if
    if (trim(restart_remap) /= 'strict' .and. trim(restart_remap) /= 'mode_index') then
      print *, "ERROR: restart_remap must be 'strict' or 'mode_index'"
      error stop 1
    end if
    if (sl_amplitude /= 0.0d0 .and. (s2_amplitude /= 0.0d0 .or. s2_mean /= 0.0d0)) then
      print *, 'ERROR: the Stokes layer and the spanwise shear S2 cannot both be active (as in hst-main)'
      error stop 1
    end if
    if (sl_amplitude /= 0.0d0 .and. sl_period <= 0.0d0) then
      print *, 'ERROR: sl_period must be positive'
      error stop 1
    end if
    if ((damp_kpar > 0.0d0) .neqv. (damp_kperp > 0.0d0)) then
      print *, 'ERROR: damp_kpar and damp_kperp must both be positive (damping on) or both zero (off)'
      error stop 1
    end if
    if (damp_norm < 1.0d0 .or. damp_edge < 0.0d0 .or. damp_edge >= 1.0d0 .or. damp_rate < 0.0d0) then
      print *, 'ERROR: damping needs damp_norm >= 1, 0 <= damp_edge < 1, damp_rate >= 0'
      error stop 1
    end if
    if (trim(initial_condition) /= 'random' .and. trim(initial_condition) /= 'taylor_green') then
      print *, "ERROR: initial_condition must be 'random' or 'taylor_green'"
      error stop 1
    end if
    if (tg_mode < 1) then
      print *, 'ERROR: tg_mode must be positive'
      error stop 1
    end if

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
    write (*, '(A,F10.6,A,F10.6,A,F8.3)') '   alfa0 =', alfa0, '   beta0 =', beta0, '   ystretch =', ystretch
    write (*, '(A,F10.2,A,F10.6,A,L1,A,L1)') '   re    =', re, '   S     =', S, '   linear =', linear, &
      '   exact_shift =', exact_shift
    if (s2_amplitude /= 0.0d0 .or. s2_mean /= 0.0d0) write (*, '(A,F10.6,A,F10.6,A,F10.6,A,F10.6)') &
      '   S2: mean=', s2_mean, '   amplitude=', s2_amplitude, '   period=', s2_period, '   start=', s2_start
    if (sl_amplitude /= 0.0d0) write (*, '(A,F8.4,A,F8.4,A,F8.4,A,F8.4,A,L1,A,L1)') &
      '   Stokes layer: amplitude=', sl_amplitude, '  period=', sl_period, '  delta=', sl_delta, &
      '  start=', sl_start, '  bodyforce=', sl_bodyforce, '  ramp=', sl_ramp
    if (damp_kpar > 0.0d0) write (*, '(A,F8.4,A,F8.4,A,F6.2,A,L1,A,F8.3,A,F6.3,A,L1,A,L1,A,F8.3)') &
      '   damping: k_par <', damp_kpar, '  k_perp <', damp_kperp, '  norm=', damp_norm, '  hard=', damp_hard, &
      '  rate=', damp_rate, '  edge=', damp_edge, '  mean=', damp_mean, '  follow=', damp_follow, '  angle=', damp_angle
    write (*, '(A,F10.6,A,F10.6,A,F10.3)') '   deltat=', deltat, '   cflmax=', cflmax, '   t_max =', t_max
    write (*, '(A,F10.4,A,F10.4,A,F10.4)') '   dt_stat=', dt_stat, '  dt_field=', dt_field, '  dt_save=', dt_save
    write (*, '(A,I10,A,F10.0,A,L1,A,L1,A,A)') '   nstep =', nstep, '   wall_max =', wall_max, &
      '   restart_from_file = ', restart_from_file, '   time_from_restart = ', &
      time_from_restart, &
       '   restart_file = ', trim(restart_file)
    write (*, '(A,A)') '   restart_remap = ', trim(restart_remap)
    write (*, '(A,L1)') '   timing = ', timing
    write (*, '(A,ES10.3,A,I6,A,F8.3)') '   amplitude=', amplitude, '   seed =', seed, '   kpeak =', kpeak
    write (*, '(A,A,A,ES10.3,A,I0)') '   initial_condition=', trim(initial_condition), &
      '   tg_amplitude=', tg_amplitude, '   tg_mode=', tg_mode
    write (*, '(A)') ' '
  end subroutine print_input

end module hst_input
