module lennard_jones
    use stdlib_kinds, only: dp
    use center_of_mass_utils, only: HardSphereType
    implicit none
    private
    public :: calculate_lj_potential

contains

  subroutine calculate_lj_potential(hard_spheres, epsilon, sigma, lj_potential, forces)
    type(HardSphereType), intent(in) :: hard_spheres
    real(dp), intent(in) :: epsilon, sigma
    real(dp), intent(out) :: lj_potential
    real(dp), intent(out) :: forces(:,:)

    integer :: i, j, num_spheres 
    real(dp) :: rij(3), r2, r6, r12, lj_pair 
    real(dp) :: force_magnitude
    real(dp), allocatable :: coords(:,:)

    num_spheres = hard_spheres%num_spheres
    coords = hard_spheres%coords

    lj_potential = 0.0_dp

    do i = 1, num_spheres - 1
        do j = i + 1, num_spheres
            rij = coords(i,:) - coords(j,:)
            r2 = dot_product(rij, rij)
            r6 = r2**3
            r12 = r6 * r6
            lj_pair = 4.0_dp * epsilon * (sigma**12 / r12 - sigma**6 / r6)
            lj_potential = lj_potential + lj_pair
                force_magnitude = -24.0_dp * epsilon * (2.0_dp * sigma**12 / r12 - sigma**6 / r6) / r2
                forces(i,:) = forces(i,:) + force_magnitude * rij
                forces(j,:) = forces(j,:) - force_magnitude * rij
        end do
    end do

  end subroutine calculate_lj_potential




end module lennard_jones
