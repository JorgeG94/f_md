module lennard_jones
   use stdlib_kinds, only: dp
   use center_of_mass_utils, only: HardSphereType
   use local_constants, only: small_distance
   use app_output
   implicit none
   private
   public :: calculate_lj_forces, calculate_lj_potential

contains

   subroutine calculate_lj_forces(hard_spheres, epsilon, sigma, lj_potential, forces)
      type(HardSphereType), intent(in) :: hard_spheres
      real(dp), intent(in) :: epsilon, sigma
      real(dp), intent(out) :: lj_potential
      real(dp), intent(out) :: forces(:,:)

      integer :: i, j, num_spheres
      real(dp) :: rij(3), r2, r6, r12, lj_pair, r
      real(dp) :: force_magnitude, forces_ij(3)
      real(dp), allocatable :: coords(:,:)
      real(dp) :: r_cutoff

      r_cutoff = 2.5_dp * sigma * 6.0_dp !* 1000.0_dp 

      num_spheres = hard_spheres%num_spheres
      coords = hard_spheres%coords


      lj_potential = 0.0_dp
      forces = 0.0_dp

      do i = 1, num_spheres - 1
         do j = i + 1, num_spheres
            rij = coords(i,:) - coords(j,:)
            r2 = dot_product(rij, rij)
            r = sqrt(r2)
            if (r < sigma) then
               print *, "distance is less than sigma", r
               print *, "Particles are too close together"
               stop
            end if
            if (r < r_cutoff) then
               lj_pair = calculate_lj_potential(rij, sigma, epsilon)
               lj_potential = lj_potential + lj_pair
               r6 = r2**3
               r12 = r6**2
               force_magnitude = 24.0_dp * epsilon * (2.0_dp * (sigma**12 / r**12) - (sigma**6 / r**6))
               forces(i,:) = forces(i,:) + force_magnitude * rij / r2
               forces(j,:) = forces(j,:) - force_magnitude * rij / r2
            end if
         end do
      end do


   end subroutine calculate_lj_forces

   function calculate_lj_potential(rij, sigma, epsilon) result(lj_pair)
      use stdlib_kinds, only: dp
      implicit none

      ! Inputs
      real(dp), intent(in) :: rij(:)      ! Distance vector between two particles
      real(dp), intent(in) :: sigma      ! Lennard-Jones sigma parameter
      real(dp), intent(in) :: epsilon    ! Lennard-Jones epsilon parameter

      ! Output
      real(dp) :: lj_pair                ! Lennard-Jones potential energy for the pair

      ! Local variables
      real(dp) :: r2, r6, r12

      ! calculate squared distance
      r2 = dot_product(rij, rij)

      ! calculate reduced distance powers
      r6 = (sigma / sqrt(r2))**6
      r12 = r6 * r6

      ! Lennard-Jones potential
      lj_pair = 4.0_dp * epsilon * (r12 - r6)
   end function calculate_lj_potential



end module lennard_jones
