!=====================================================================!
!     postpro: plane statistics of the snapshots of a run             !
!=====================================================================!
!
! Usage:  mpirun -np N postpro [postpro.in]
!
! postpro.in is a namelist (the copy in the repository is commented) that
! names the snapshots, the deck of the run and what to compute:
!
!   &postpro  first = 1, last = 20, step = 1,  deck = 'hst.in',  outdir = 'statistics',
!             mean = .true., stresses = .true., spectra = .true.,  budgets = 'uu vv ww uv uw vw' /
!
! Reads fields/field<n>.fld for n = first, first+step, ..., last and writes
! into outdir the averages over these snapshots as functions of y, one line
! per y point, plain text with a header line:
!
!   mean.dat        y  U  W  P        the plane means (the (0,0) mode; the mean
!                                     shear S*y itself is not included)
!   stresses.dat    y  uu vv ww uv uw vw        the Reynolds stresses <u_i u_j>
!   budget_uu.dat   y  prod diss pstrain ttrsp ptrsp vdiff sum
!   ... one file per component named in `budgets` (any of uu vv ww uv uw vw;
!         an empty string skips the pressure and the gradients, the expensive
!         part): the terms of
!         d<u_i u_j>/dt = prod + pstrain - diss + ttrsp + ptrsp + vdiff  (= sum,
!                         zero in a statistically stationary flow)
!         prod    = -<u_i v> dU_j/dy - <u_j v> dU_i/dy  (U = S y + U(y), W = S2 y + W(y))
!         diss    = 2 nu <du_i/dx_k du_j/dx_k>
!         pstrain = <p (du_i/dx_j + du_j/dx_i)>
!         ttrsp   = -d/dy <v u_i u_j>
!         ptrsp   = -d/dy <p u_i delta_jv + p u_j delta_iv>
!         vdiff   = nu d2<u_i u_j>/dy2
!   spectra_xz.dat  kx kz Euu Evv Eww Euv Euw Evw    two-dimensional spectra,
!                   averaged over y, -kz folded onto kz (summed over kx, kz
!                   they give the y-averaged stresses)
!   spectra_x.dat   kx Euu ... (summed over kz);  spectra_z.dat  kz Euu ... (over kx)
!
! Axes as in the solver: x streamwise, y the shear direction (velocity v),
! z spanwise (w).  Fluctuations are taken about the plane mean of each
! snapshot.  The pressure is recomputed from the velocity (hst_pressure),
! so p_fields/ is not needed.  Everything else is the solver's: the file
! reader, the shear-periodic ghost rows, the compact d/dy, the transforms to
! physical space.  The plane averages of products are sums over the
! dealiased physical grid (exact for a product of two fields, aliased by
! the 3/2 padding for the triple products, as in the CPL tools), and the
! y derivatives of the averaged profiles use the solver's compact stencils
! (D0 f' = D1 f, a small dense periodic system solved once and for all).
program hst_postpro

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
  use hst_io
  use hst_pressure

  implicit none

  ! the fields taken to physical space: u v w p, then du_i/dx_j (i outer)
  integer, parameter :: NF = 13
  integer, parameter :: IU = 1, IV = 2, IW = 3, IP = 4
  ! the six stress components ij and their velocity indices
  integer, parameter :: CI(6) = [1, 2, 3, 1, 1, 2], CJ(6) = [1, 2, 3, 2, 3, 3]
  character(len=2), parameter :: NAME(6) = ['uu', 'vv', 'ww', 'uv', 'uw', 'vw']

  ! the input file
  integer :: first, last, step
  character(len=256) :: deck, outdir, budgets
  logical :: mean, stresses, spectra
  namelist /postpro/ first, last, step, deck, outdir, mean, stresses, spectra, budgets

  integer :: n, nsnap, ierr, c, iy, unit, ios
  character(len=256) :: arg
  character(len=64) :: fname
  logical :: exists, want_budget(6), any_budget

  ! host copies with the layout of V: the fluctuations, the pressure, d/dy
  complex(C_DOUBLE_COMPLEX), allocatable :: Vf(:, :, :, :), p(:, :, :), dVdy(:, :, :, :)
  ! the NF fields in physical space, this rank's x points, z lines and y rows
  real(C_DOUBLE), allocatable :: phys(:, :, :, :)
  ! accumulated profiles (all y rows; a rank adds its own rows, the sum over ranks completes them)
  real(C_DOUBLE), allocatable :: umean_sum(:, :), stress(:, :), prod(:, :), diss(:, :), pstrain(:, :), &
                                 tflux(:, :), pflux(:, :), spec(:, :, :)
  real(C_DOUBLE), allocatable :: umean(:, :), dUdy(:, :)     ! this snapshot's mean profiles and their slopes
  real(C_DOUBLE), allocatable :: D0lu(:, :)                  ! the compact D0 of the y profiles, LU-factorized

  !------------------------------------------------------------ set-up ----
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, iproc, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
  has_terminal = (iproc == 0)
  call select_device()
  ! the input file: defaults, then the namelist
  first = 1; last = 1; step = 1; deck = 'hst.in'; outdir = 'statistics'
  mean = .true.; stresses = .true.; spectra = .true.; budgets = 'uu vv ww uv uw vw'
  arg = 'postpro.in'
  if (command_argument_count() >= 1) call get_command_argument(1, arg)
  open (newunit=unit, file=trim(arg), status='old', action='read', iostat=ios)
  if (ios /= 0) then
    if (has_terminal) print *, 'ERROR: cannot open '//trim(arg)//'   (usage: mpirun -np N postpro [postpro.in])'
    error stop 1
  end if
  read (unit, nml=postpro)
  close (unit)
  want_budget = [(index(budgets, NAME(c)) > 0, c=1, 6)]
  any_budget = any(want_budget)
  call read_input(trim(deck))
  time_from_restart = .true.           ! the clock of each snapshot comes from its file
  call setup_decomposition()
  call allocate_fields()
  call init_fft()
  call init_linsolve()
  call setup_derivatives()
  if (has_terminal) print '(A,I0,A,I0,A,I0,A,I0,A,A,A,L1,A,L1,A,L1,A,A)', ' postpro: snapshots ', first, '..', last, &
    ' by ', step, ' on ', nproc, ' ranks; deck ', trim(deck), '; mean ', mean, ', stresses ', stresses, &
    ', spectra ', spectra, ', budgets ', trim(budgets)

  allocate (Vf, mold=V)
  allocate (p(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN), dVdy(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN, 3))
  allocate (phys(2*nxd, nzB, ny0:nyN, NF))
  allocate (umean_sum(0:ny - 1, 3), stress(0:ny - 1, 6), prod(0:ny - 1, 6), diss(0:ny - 1, 6), &
            pstrain(0:ny - 1, 6), tflux(0:ny - 1, 6), pflux(0:ny - 1, 6), spec(0:nx, 0:nz, 6))
  allocate (umean(0:ny - 1, 3), dUdy(0:ny - 1, 3))
  call setup_profile_derivatives()
  umean_sum = 0; stress = 0; prod = 0; diss = 0; pstrain = 0; tflux = 0; pflux = 0; spec = 0
  nsnap = 0

  !------------------------------------------------ loop over snapshots ----
  do n = first, last, step
    write (fname, '(A,I0,A)') 'fields/field', n, '.fld'
    inquire (file=fname, exist=exists)
    if (.not. exists) then
      if (has_terminal) print *, ' ', trim(fname), ' not found: stopping here'
      exit
    end if
    call restart_read(fname)                     ! V on the host, time from the file
    !$omp target update to(V)
    do c = 1, 3
      call fill_ghosts(c)
    end do
    nsnap = nsnap + 1
    if (has_terminal) print '(A,A,A,F10.4)', '   ', trim(fname), ' at time', time

    ! the budgets need the pressure (device, rhs(:, :, :, 2)) and d/dy of the
    ! velocity with the compact scheme through the ghost rows
    p = 0; dVdy = 0
    if (any_budget) then
      call compute_pressure(rhs(:, :, :, 2))
      !$omp target update from(rhs)
      p = rhs(:, :, :, 2)
      do c = 1, 3
        call line_solve(KIND_DY, 0.0d0, V(:, :, :, c), rhs(:, :, :, 1))
        !$omp target update from(rhs)
        dVdy(:, :, :, c) = rhs(:, :, :, 1)
      end do
    end if

    ! plane means (the (0,0) mode, real) of this snapshot and their slopes;
    ! the fluctuations are everything else
    umean = 0
    if (nx0 == 0) then
      do iy = ny0, nyN
        umean(iy, :) = [dreal(V(iy, 0, 0, 1)), dreal(V(iy, 0, 0, 3)), dreal(p(iy, 0, 0))]
      end do
    end if
    call MPI_Allreduce(MPI_IN_PLACE, umean, size(umean), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    umean_sum = umean_sum + umean
    dUdy(:, 1) = S + ddy(umean(:, 1))            ! dU/dy: the mean shear plus the profile's
    dUdy(:, 2) = 0.0d0
    dUdy(:, 3) = s2_of(time) + ddy(umean(:, 2))  ! dW/dy
    Vf = V
    if (nx0 == 0) then
      Vf(:, 0, 0, :) = 0; p(:, 0, 0) = 0; dVdy(:, 0, 0, :) = 0
    end if

    call accumulate_spectral()
    if (any_budget) then
      call to_physical_space()
      call accumulate_physical()
    end if
  end do

  !--------------------------------------------------- average and write ----
  if (nsnap == 0) then
    if (has_terminal) print *, ' no snapshot read'
    call MPI_Finalize(ierr); stop
  end if
  call MPI_Allreduce(MPI_IN_PLACE, stress, size(stress), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_Allreduce(MPI_IN_PLACE, prod, size(prod), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_Allreduce(MPI_IN_PLACE, diss, size(diss), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_Allreduce(MPI_IN_PLACE, pstrain, size(pstrain), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_Allreduce(MPI_IN_PLACE, tflux, size(tflux), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_Allreduce(MPI_IN_PLACE, pflux, size(pflux), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_Allreduce(MPI_IN_PLACE, spec, size(spec), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  umean_sum = umean_sum/nsnap; stress = stress/nsnap; prod = prod/nsnap; spec = spec/nsnap
  ! the physical-space sums: per snapshot and per point of the dealiased plane
  diss = diss/(nsnap*2.0d0*nxd*nzd); pstrain = pstrain/(nsnap*2.0d0*nxd*nzd)
  tflux = tflux/(nsnap*2.0d0*nxd*nzd); pflux = pflux/(nsnap*2.0d0*nxd*nzd)
  if (has_terminal) call write_statistics()

  call free_linsolve()
  call free_fft()
  call free_fields()
  call free_mpi()
  call MPI_Finalize(ierr)

contains

  ! Stresses per plane, their production, and the spectra: sums over this
  ! rank's modes (ix > 0 stands for -ix too) and rows.
  subroutine accumulate_spectral()
    integer :: ix, iz, iy, ij, a, b
    real(C_DOUBLE) :: w, s
    do ix = nx0, nxN
      do iz = -nz, nz
        if (ix == 0 .and. iz == 0) cycle
        w = 2.0d0; if (ix == 0) w = 1.0d0
        do iy = ny0, nyN
          do ij = 1, 6
            a = CI(ij); b = CJ(ij)
            s = w*dreal(Vf(iy, iz, ix, a)*conjg(Vf(iy, iz, ix, b)))
            stress(iy, ij) = stress(iy, ij) + s
            spec(ix, abs(iz), ij) = spec(ix, abs(iz), ij) + s*dyl(iy)/ly
            prod(iy, ij) = prod(iy, ij) - w*dreal(Vf(iy, iz, ix, a)*conjg(Vf(iy, iz, ix, IV)))*dUdy(iy, b) &
                                        - w*dreal(Vf(iy, iz, ix, b)*conjg(Vf(iy, iz, ix, IV)))*dUdy(iy, a)
          end do
        end do
      end do
    end do
  end subroutine accumulate_spectral

  ! The NF fields to physical space, three at a time through V and the
  ! solver's transform (its output rVVdx is the physical field on the
  ! dealiased grid, this rank's z lines).
  subroutine to_physical_space()
    integer :: g, m, f
    do g = 1, (NF + 2)/3
      do m = 1, 3
        f = 3*(g - 1) + m
        if (f <= NF) then
          call spectral_field(f, V(:, :, :, m))
        else
          V(:, :, :, m) = 0
        end if
      end do
      !$omp target update to(V)
      call transform_to_physical()
      !$omp target update from(VVdx)
      do m = 1, 3
        f = 3*(g - 1) + m
        if (f <= NF) phys(:, :, :, f) = rVVdx(1:2*nxd, :, ny0:nyN, m)
      end do
    end do
  end subroutine to_physical_space

  ! Field f in spectral space: u, v, w, p, or du_a/dx_b with x, z by the
  ! wavenumbers and y from the compact derivative.
  subroutine spectral_field(f, out)
    integer, intent(in) :: f
    complex(C_DOUBLE_COMPLEX), intent(out) :: out(ny0 - 2:, -nz:, nx0:)
    integer :: a, b, ix, iz
    if (f <= 3) then
      out = Vf(:, :, :, f)
    else if (f == IP) then
      out = p
    else
      a = (f - IP - 1)/3 + 1; b = mod(f - IP - 1, 3) + 1
      select case (b)
      case (1)
        do ix = nx0, nxN
          out(:, :, ix) = ialfa(ix)*Vf(:, :, ix, a)
        end do
      case (2)
        out = dVdy(:, :, :, a)
      case (3)
        do iz = -nz, nz
          out(:, iz, :) = ibeta(iz)*Vf(:, iz, :, a)
        end do
      end select
    end if
  end subroutine spectral_field

  ! Plane sums of the products that need physical space: the dissipation
  ! tensor, the pressure-strain, the triple-product and pressure fluxes.
  subroutine accumulate_physical()
    integer :: i, k, iy, ij, a, b, l
    real(C_DOUBLE) :: u(3), pp, grad(3, 3)
    do iy = ny0, nyN
      do k = 1, nzB
        do i = 1, 2*nxd
          u = phys(i, k, iy, IU:IW); pp = phys(i, k, iy, IP)
          grad = reshape(phys(i, k, iy, IP + 1:NF), [3, 3], order=[2, 1])   ! grad(a, b) = du_a/dx_b
          do ij = 1, 6
            a = CI(ij); b = CJ(ij)
            diss(iy, ij) = diss(iy, ij) + 2.0d0*ni*sum(grad(a, :)*grad(b, :))
            pstrain(iy, ij) = pstrain(iy, ij) + pp*(grad(a, b) + grad(b, a))
            tflux(iy, ij) = tflux(iy, ij) + u(IV)*u(a)*u(b)
            if (b == IV) pflux(iy, ij) = pflux(iy, ij) + pp*u(a)
            if (a == IV) pflux(iy, ij) = pflux(iy, ij) + pp*u(b)
          end do
        end do
      end do
    end do
  end subroutine accumulate_physical

  ! d/dy and d2/dy2 of a periodic profile with the solver's compact
  ! stencils der(iy, k, -2:2) (hst_derivatives): D0 f' = D1 f and D0 f'' =
  ! D2 f, where D0 is a periodic pentadiagonal matrix.  It is small (ny x
  ! ny), so it is stored dense and LU-factorized once, without pivoting
  ! (it is diagonally dominant).
  subroutine setup_profile_derivatives()
    integer :: iy, j, k
    allocate (D0lu(0:ny - 1, 0:ny - 1)); D0lu = 0
    do iy = 0, ny - 1
      do j = -2, 2
        D0lu(iy, modulo(iy + j, ny)) = D0lu(iy, modulo(iy + j, ny)) + der(iy, 0, j)
      end do
    end do
    do k = 0, ny - 2                       ! Gaussian elimination, the multipliers kept below the diagonal
      do iy = k + 1, ny - 1
        D0lu(iy, k) = D0lu(iy, k)/D0lu(k, k)
        D0lu(iy, k + 1:) = D0lu(iy, k + 1:) - D0lu(iy, k)*D0lu(k, k + 1:)
      end do
    end do
  end subroutine setup_profile_derivatives

  ! The compact derivative of order k (1 or 2) of the periodic profile f.
  function compact_derivative(f, k) result(d)
    real(C_DOUBLE), intent(in) :: f(0:ny - 1)
    integer, intent(in) :: k
    real(C_DOUBLE) :: d(0:ny - 1)
    integer :: iy, j
    do iy = 0, ny - 1                      ! the right-hand side Dk f
      d(iy) = sum([(der(iy, k, j)*f(modulo(iy + j, ny)), j=-2, 2)])
    end do
    do iy = 1, ny - 1                      ! forward substitution, then back substitution
      d(iy) = d(iy) - sum(D0lu(iy, 0:iy - 1)*d(0:iy - 1))
    end do
    do iy = ny - 1, 0, -1
      d(iy) = (d(iy) - sum(D0lu(iy, iy + 1:)*d(iy + 1:)))/D0lu(iy, iy)
    end do
  end function compact_derivative

  function ddy(f) result(d)
    real(C_DOUBLE), intent(in) :: f(0:ny - 1)
    real(C_DOUBLE) :: d(0:ny - 1)
    d = compact_derivative(f, 1)
  end function ddy

  function d2dy2(f) result(d)
    real(C_DOUBLE), intent(in) :: f(0:ny - 1)
    real(C_DOUBLE) :: d(0:ny - 1)
    d = compact_derivative(f, 2)
  end function d2dy2

  ! The files of statistics/ and a summary on the screen.
  subroutine write_statistics()
    integer :: unit, iy, ij, ix, iz
    real(C_DOUBLE) :: ttrsp(0:ny - 1), ptrsp(0:ny - 1), vdiff(0:ny - 1), total(0:ny - 1)
    real(C_DOUBLE) :: q2, box(6, 7)
    character(len=*), parameter :: fmt = '(*(ES16.8, 1X))'

    call execute_command_line('mkdir -p '//trim(outdir))
    if (mean) then
      open (newunit=unit, file=trim(outdir)//'/mean.dat', action='write')
      write (unit, '(A)') '# y  U  W  P   (plane means; the mean flow S y is not included)'
      do iy = 0, ny - 1
        write (unit, fmt) y(iy), umean_sum(iy, :)
      end do
      close (unit)
    end if

    if (stresses) then
      open (newunit=unit, file=trim(outdir)//'/stresses.dat', action='write')
      write (unit, '(A)') '# y  uu  vv  ww  uv  uw  vw'
      do iy = 0, ny - 1
        write (unit, fmt) y(iy), stress(iy, :)
      end do
      close (unit)
    end if

    box = 0
    do ij = 1, 6
      box(ij, 1) = sum(stress(:, ij)*dyl)/ly
      if (.not. want_budget(ij)) cycle
      ttrsp = -ddy(tflux(:, ij)); ptrsp = -ddy(pflux(:, ij)); vdiff = ni*d2dy2(stress(:, ij))
      total = prod(:, ij) + pstrain(:, ij) - diss(:, ij) + ttrsp + ptrsp + vdiff
      open (newunit=unit, file=trim(outdir)//'/budget_'//NAME(ij)//'.dat', action='write')
      write (unit, '(A)') '# y  prod  diss  pstrain  ttrsp  ptrsp  vdiff  sum      (d<'//NAME(ij)// &
        '>/dt = prod + pstrain - diss + ttrsp + ptrsp + vdiff)'
      do iy = 0, ny - 1
        write (unit, fmt) y(iy), prod(iy, ij), diss(iy, ij), pstrain(iy, ij), ttrsp(iy), ptrsp(iy), vdiff(iy), total(iy)
      end do
      close (unit)
      ! box averages for the summary
      box(ij, 2:7) = [sum(prod(:, ij)*dyl), sum(diss(:, ij)*dyl), sum(pstrain(:, ij)*dyl), &
                      sum(ttrsp*dyl), sum(ptrsp*dyl), sum(vdiff*dyl)]/ly
    end do

    if (spectra) then
    open (newunit=unit, file=trim(outdir)//'/spectra_xz.dat', action='write')
    write (unit, '(A)') '# kx  kz  Euu  Evv  Eww  Euv  Euw  Evw   (y-averaged, kz >= 0; the sum over kx, kz is <u_i u_j>)'
    do ix = 0, nx
      do iz = 0, nz
        write (unit, fmt) alfa0*ix, beta0*iz, spec(ix, iz, :)
      end do
    end do
    close (unit)
    open (newunit=unit, file=trim(outdir)//'/spectra_x.dat', action='write')
    write (unit, '(A)') '# kx  Euu  Evv  Eww  Euv  Euw  Evw   (summed over kz)'
    do ix = 0, nx
      write (unit, fmt) alfa0*ix, sum(spec(ix, :, :), dim=1)
    end do
    close (unit)
    open (newunit=unit, file=trim(outdir)//'/spectra_z.dat', action='write')
    write (unit, '(A)') '# kz  Euu  Evv  Eww  Euv  Euw  Evw   (summed over kx)'
    do iz = 0, nz
      write (unit, fmt) beta0*iz, sum(spec(:, iz, :), dim=1)
    end do
    close (unit)
    end if

    q2 = box(1, 1) + box(2, 1) + box(3, 1)
    print '(A,A,A,I0,A)', ' ', trim(outdir), '/ written: averages over ', nsnap, ' snapshots'
    print '(A,F10.5,A,F8.4)', '   q2 =', q2, '   -uv/q2 =', -box(4, 1)/q2
    if (all(want_budget(1:3))) print '(A,F10.5,A,F8.3)', '   eps =', 0.5d0*(box(1, 3) + box(2, 3) + box(3, 3)), &
      '   S* =', S*q2/(0.5d0*(box(1, 3) + box(2, 3) + box(3, 3)))
    if (any_budget) print '(A)', &
      '   box budget:   <ij>        prod        diss     pstrain       ttrsp       ptrsp       vdiff         sum'
    do ij = 1, 6
      if (want_budget(ij)) print '(A,A,8(1X,ES11.3))', '     ', NAME(ij), box(ij, 1), box(ij, 2), box(ij, 3), box(ij, 4), &
        box(ij, 5), box(ij, 6), box(ij, 7), box(ij, 2) + box(ij, 4) - box(ij, 3) + box(ij, 5) + box(ij, 6) + box(ij, 7)
    end do
  end subroutine write_statistics

end program hst_postpro
