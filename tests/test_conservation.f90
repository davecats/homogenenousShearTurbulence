! Energy conservation of the nonlinear terms, product by product (WP5).
!
!   mpirun -np N build-cpu/test_conservation [hst.in]     (S = 0, inviscid, forced)
!
! One inviscid step of length dt from the seeded random divergence-free
! field, with only one of the six products (uu, vv, ww, uv, vw, uw) kept,
! and the energy change per unit time divided by the energy.  In the
! continuum each of uu, vv, ww alone injects nothing (u d_x(u^2) integrates
! to zero), the three cross products sum to zero, and so does the whole.
! Discretely the numbers are O(dt) plus the (small) non-conservation of the
! compact scheme; anything of order one marks a wrong term.  Use a deck with
! a well-resolved field (tests/decks/pressure.in).
program test_conservation

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
  use hst_initial
#ifdef HAVE_CUDA
  use omp_lib
#endif

  implicit none

  character(len=256) :: deck
  integer :: ierr, mask, i, m, j, kk, ii, y_first, y_last
  real(C_DOUBLE) :: q0, q1, rate(0:6), dt_sub
  character(len=2), parameter :: name(6) = ['uu', 'vv', 'ww', 'uv', 'vw', 'uw']

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, iproc, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
  has_terminal = (iproc == 0)
#ifdef HAVE_CUDA
  call omp_set_default_device(mod(iproc, omp_get_num_devices()))
#endif
  deck = 'hst.in'
  if (command_argument_count() >= 1) call get_command_argument(1, deck)
  call read_input(trim(deck))
  S = 0.0d0; linear = .false.; ni = 1.0d-12
  call setup_decomposition()
  call allocate_fields()
  call init_fft()
  call init_linsolve()
  call setup_derivatives()
  y_first = ny0 - 2
  y_last = nyN + 2

  do mask = 0, 6
    time = 0.0d0; deltat = 1.0d-4; oldrhs = 0
    !$omp target update to(oldrhs)
    call generate_initial_field()
    !$omp target update to(V)
    do m = 1, 3
      call fill_ghosts(m)
    end do
    ! one unmasked step first: the solver projects the analytic field onto
    ! its discretely divergence-free space, which changes the energy at
    ! O(dy^6) -- large against dt.  Measure from the projected field.
    call timestep()
    q0 = energy()
    ! one RK step, as timestep(), with the product mask
    do i = 1, 3
      call transform_to_physical()
      call buildrhs_prepare(RK_rai(:, i))
      do m = 1, 6
        call build_products(m)
        if (mask /= 0 .and. m /= mask) then
          !$omp target teams distribute parallel do collapse(3) default(none) &
          !$omp shared(products, nxd, nzB, y_first, y_last) private(j, kk, ii)
          do j = y_first, y_last
            do kk = 1, nzB
              do ii = 1, 2*nxd
                products(ii, kk, j) = 0.0d0
              end do
            end do
          end do
        end if
        call products_to_spectral()
        call buildrhs(RK_rai(:, i), m)
      end do
      dt_sub = 2.0d0/RK_rai(1, i)*deltat
      call shear_shift(dt_sub)
      time = time + dt_sub
      call linsolve(RK_rai(1, i)/deltat)
    end do
    q1 = energy()
    rate(mask) = (q1 - q0)/deltat/q0
    if (has_terminal) then
      if (mask == 0) then
        write (*, '(A,ES11.3)') '   all products : (dq2/dt)/q2 = ', rate(mask)
      else
        write (*, '(A,A,A,ES11.3)') '   ', name(mask), ' alone     : (dq2/dt)/q2 = ', rate(mask)
      end if
    end if
  end do
  if (has_terminal) write (*, '(A,ES11.3)') '   uv + vw + uw : (dq2/dt)/q2 = ', rate(4) + rate(5) + rate(6)

  call free_linsolve(); call free_fft(); call free_fields(); call free_mpi()
  call MPI_Finalize(ierr)
  if (maxval(abs(rate(0:3))) > 1.0d-2 .or. abs(rate(4) + rate(5) + rate(6)) > 1.0d-2) then
    print *, 'FAILED'
    error stop 1
  end if
  if (has_terminal) print *, 'PASSED'

contains

  ! q2 = <u_i u_i>, all ranks
  real(C_DOUBLE) function energy()
    integer :: ix, iz, iy, c
    real(C_DOUBLE) :: s, w, g
    !$omp target update from(V)
    s = 0.0d0
    do ix = nx0, nxN
      w = 2.0d0
      if (ix == 0) w = 1.0d0
      do iz = -nz, nz
        do iy = 0, ny - 1
          do c = 1, 3
            s = s + w*abs(V(iy, iz, ix, c))**2
          end do
        end do
      end do
    end do
    call MPI_Allreduce(s, g, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    energy = g/ny
  end function energy

end program test_conservation
