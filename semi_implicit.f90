program semi_implicit
  implicit none

  integer, parameter :: dp = kind(0.0d0)

  real(dp), parameter :: e = 1.602176634e-19_dp
  real(dp), parameter :: eps0 = 8.8541878176e-12_dp

  integer, parameter  :: nx        = 300
  real(dp), parameter :: L         = 10e-3_dp
  real(dp), parameter :: dx        = L / nx
  real(dp), parameter :: t_end     = 30e-9_dp
  real(dp), parameter :: dt_output = 1e-9_dp
  real(dp), parameter :: mu        = 0.03_dp
  real(dp), parameter :: D         = 0.1_dp
  real(dp), parameter :: E0        = 1e6_dp
  real(dp), parameter :: n0        = 1e20_dp

  real(dp) :: dt, dt_cfl, dt_diff, dt_tau, dt_max

  real(dp) :: x(-1:nx+2)
  real(dp) :: n_e(-1:nx+2)
  real(dp) :: n_i(-1:nx+2)
  real(dp) :: field_fc(nx+1)
  real(dp) :: field_cc(nx)
  real(dp) :: phi(nx)
  real(dp) :: time
  integer  :: n, n_output
  logical  :: output_this_step

  dt_cfl  = 0.5_dp * dx / max(abs(mu * E0), 1e-10_dp)
  dt_diff = 0.5_dp * dx**2 / max(D, 1e-10_dp)
  dt_tau  = eps0 / (e * mu * n0)
  dt_max  = 1/(1/dt_cfl + 1/dt_diff)

  print *, "L:         ", L
  print *, "nx:        ", nx
  print *, "dx:        ", dx
  print *, "t_end:     ", t_end
  print *, "dt_output: ", dt_output
  print *, "mu:        ", mu
  print *, "D:         ", D
  print *, "E0:        ", E0
  print *, "dt:        ", dt
  print *, "tau:       ", dt_tau

  do n = -1, nx+2
     x(n) = (n - 0.5_dp) * dx
  end do

  ! Initial conditions
  where (x > 0.4_dp * L .and. x < 0.6 * L)
     n_e = n0
     n_i = n0
  elsewhere
     n_e = 0
     n_i = 0
  end where

  call solve_field(0.0_dp)

  time     = 0.0_dp
  n_output = 0

  do while (time <= t_end)
     output_this_step = (time + dt_max >= dt_output * n_output)

     if (output_this_step) then
        dt = dt_output * n_output - time
        n_output = n_output + 1
     else
        dt = dt_max
     end if

     call solve_field(dt)
     call advance(dt)
     time = time + dt

     if (output_this_step) then
        call write_output(n_output)
     end if
  end do

contains

  subroutine solve_field(dt)
    real(dp), intent(in)    :: dt
    real(dp)                :: n_face(nx+1)
    real(dp)                :: coeff(nx+1)
    real(dp), dimension(nx) :: l_diag, diag, u_diag, rhs
    real(dp)                :: phi_left, phi_right
    integer                 :: i, n

    do n = 1, nx + 1
       ! Determine n_e at cell faces
       n_face(n) = 0.5_dp * (n_e(n-1) + n_e(n))
       ! Set effective epsilon for Poisson equation
       coeff(n) = -eps0 / e - dt * mu * n_face(n)
    end do

    ! TODO: maybe check that densities are zero at boundaries?
    l_diag(1:nx) = coeff(1:nx) / dx**2
    ! For the left boundary condition at the cell face
    l_diag(1) = 2 * l_diag(1)
    u_diag(1:nx) = coeff(2:nx+1) / dx**2
    ! For the right boundary condition at the cell face
    u_diag(nx) = 2 * u_diag(nx)
    diag = -(l_diag+u_diag)

    ! Right-hand side
    rhs = n_i(1:nx) - n_e(1:nx)

    ! Add diffusive terms to rhs
    do i = 1, nx
       rhs(i) = rhs(i) - dt * D/dx**2 * (n_e(i-1) - 2 * n_e(i) + n_e(i+1))
    end do

    ! Set Dirichlet boundary conditions
    phi_left = 0.0_dp
    phi_right = -E0 * L
    rhs(1) = rhs(1) - phi_left * l_diag(1)
    rhs(nx) = rhs(nx) - phi_right * u_diag(nx)

    call solve_tridiag(l_diag, diag, u_diag, rhs, phi, nx)

    field_fc(1) = -(phi(1) - phi_left) / (0.5_dp * dx)
    field_fc(nx+1) = -(phi_right - phi(nx)) / (0.5_dp * dx)
    do i = 2, nx
       field_fc(i) = -(phi(i) - phi(i-1)) / dx
    end do

    field_cc = 0.5_dp * (field_fc(2:) + field_fc(1:nx))

  end subroutine solve_field

  subroutine advance(dt)
    real(dp), intent(in) :: dt
    real(dp)             :: flux(1:nx+1), Dvec(1:nx+1)
    integer              :: i

    Dvec = D
    call get_flux_1d(n_e, -field_fc * mu, Dvec, dx, flux, nx, 2)

    do i = 1, nx
       n_e(i) = n_e(i) + dt * (flux(i) - flux(i+1))/dx
    end do
  end subroutine advance

  ! Solve a tridiagonal matrix
  ! Copied from:
  ! https://en.wikibooks.org/wiki/Algorithm_Implementation/
  ! Linear_Algebra/Tridiagonal_matrix_algorithm#Fortran_90
  subroutine solve_tridiag(a, b, c, d, x, n)
    !> number of equations
    integer, intent(in)    :: n
    !> sub-diagonal, a(2:n) is used
    real(dp), intent(in)   :: a(n)
    !> the main diagonal
    real(dp), intent(in)   :: b(n)
    !> sup-diagonal, c(1:n-1) is used
    real(dp), intent(in)   :: c(n)
    !> right-hand side
    real(dp), intent(in)   :: d(n)
    !> solution
    real(dp), intent(out)  :: x(n)
    real(dp), dimension(n) :: cc, dd
    real(dp)               :: m
    integer                :: i

    ! initialize c-prime and d-prime
    cc(1) = c(1)/b(1)
    dd(1) = d(1)/b(1)
    ! solve for vectors c-prime and d-prime
    do i = 2,n
       m = b(i)-cc(i-1)*a(i)
       cc(i) = c(i)/m
       dd(i) = (d(i)-dd(i-1)*a(i))/m
    enddo
    ! initialize x
    x(n) = dd(n)
    ! solve for x from the vectors c-prime and d-prime
    do i = n-1, 1, -1
       x(i) = dd(i)-cc(i)*x(i+1)
    end do

  end subroutine solve_tridiag

  !> Compute advective and diffusive flux
  subroutine get_flux_1d(cc, v, dc, dx, flux, nc, ngc)
    integer, intent(in)   :: nc               !< Number of cells
    integer, intent(in)   :: ngc              !< Number of ghost cells
    real(dp), intent(in)  :: cc(1-ngc:nc+ngc) !< Cell-centered values
    !> Input: velocities at cell faces
    real(dp), intent(in)  :: v(1:nc+1)
    !> Input: diffusion coefficients at cell faces
    real(dp), intent(in)  :: dc(1:nc+1)
    !> Grid spacing
    real(dp), intent(in)  :: dx
    !> Output: flux at cell faces
    real(dp), intent(out) :: flux(1:nc+1)
    real(dp)              :: gradp, gradc, gradn, inv_dx
    integer               :: n

    inv_dx = 1/dx

    do n = 1, nc+1
       gradc = cc(n) - cc(n-1)  ! Current gradient
       if (v(n) < 0.0_dp) then
          gradn = cc(n+1) - cc(n) ! Next gradient
          flux(n) = v(n) * (cc(n) - koren_mlim(gradc, gradn))
       else                     ! v(n) > 0
          gradp = cc(n-1) - cc(n-2) ! Previous gradient
          flux(n) = v(n) * (cc(n-1) + koren_mlim(gradc, gradp))
       end if
       ! Add diffusive flux (central differences)
       flux(n) = flux(n) - dc(n) * gradc * inv_dx
    end do

  end subroutine get_flux_1d

  !> Modified implementation of Koren limiter, to avoid division and the min/max
  !> functions, which can be problematic / expensive. In most literature, you
  !> have r = a / b (ratio of gradients). Then the limiter phi(r) is multiplied
  !> with b. With this implementation, you get phi(r) * b
  elemental function koren_mlim(a, b) result(bphi)
    real(dp), intent(in) :: a  !< Density gradient (numerator)
    real(dp), intent(in) :: b  !< Density gradient (denominator)
    real(dp), parameter  :: sixth = 1/6.0_dp
    real(dp)             :: bphi, aa, ab

    aa = a * a
    ab = a * b

    if (ab <= 0) then
       ! a and b have different sign or one of them is zero, so r is either 0,
       ! inf or negative (special case a == b == 0 is ignored)
       bphi = 0
    else if (aa <= 0.25_dp * ab) then
       ! 0 < a/b <= 1/4, limiter has value a/b
       bphi = a
    else if (aa <= 2.5_dp * ab) then
       ! 1/4 < a/b <= 2.5, limiter has value (1+2*a/b)/6
       bphi = sixth * (b + 2*a)
    else
       ! (1+2*a/b)/6 >= 1, limiter has value 1
       bphi = b
    end if
  end function koren_mlim

  subroutine write_output(ix)
    integer, intent(in) :: ix
    character(len=200)  :: fname
    integer             :: n, my_unit

    write(fname, "(A,I0.6,A)") "output_", ix, ".txt"
    open(newunit=my_unit, file=trim(fname))

    write(my_unit, "(A)") "x field electron pos_ion potential"
    do n = 1, nx
       write(my_unit, *) x(n), &
            field_cc(n), &
            n_e(n), &
            n_i(n), &
            phi(n)
    end do
    close(my_unit)

    print *, "Written ", trim(fname)
  end subroutine write_output

end program
