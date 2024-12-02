program check
  use types_module
  use lennard_jones
implicit none    ! Parameters
    real(dp) :: sigma, epsilon, r_min, r_max, r_step, r, lj_potential
    integer :: num_points, i

    ! File for output
    character(len=100) :: output_file
    output_file = "lj_potential.dat"

    ! Input parameters
    sigma = 1.0_dp       ! Reduced sigma
    epsilon = 1.0_dp     ! Reduced epsilon
    r_min = 0.8_dp       ! Start distance (must be > 0 to avoid divide-by-zero)
    r_max = 3.0_dp       ! End distance
    r_step = 0.01_dp     ! Step size

    ! Number of points
    num_points = int((r_max - r_min) / r_step) + 1

    ! Open output file
    open(unit=10, file=output_file, status='replace')

    ! Compute and write Lennard-Jones potential
    write(10, '(A)') "# r, LJ Potential"
    do i = 0, num_points - 1
        r = r_min + i * r_step
        lj_potential = calculate_lj_potential([r, 0.0_dp, 0.0_dp], sigma, epsilon)
        write(10, '(F10.4, 1X, F12.6)') r, lj_potential
    end do

    close(10)
    print *, "LJ potential saved to:", trim(output_file)

print *, "Put some tests in here!"
end program check
