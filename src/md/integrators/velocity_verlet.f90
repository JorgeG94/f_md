module velocity_verlet
    use stdlib_kinds, only: dp
    use center_of_mass_utils, only: HardSphereType
    implicit none
    private 
    public :: do_velocity_verlet_step 
    contains

    subroutine do_velocity_verlet_step(hard_spheres, forces, dt, masses, velocities)
        type(HardSphereType), intent(inout) :: hard_spheres
        real(dp), intent(in) :: forces(:,:)
        real(dp), intent(in) :: dt
        real(dp), intent(in) :: masses(:)
        real(dp), intent(inout) :: velocities(:,:)
        
        integer :: i, num_spheres

        num_spheres = hard_spheres%num_spheres

        do i = 1, num_spheres 
          hard_spheres%coords(i,:) = hard_spheres%coords(i,:) + &
          velocities(i,:) * dt + 0.5_dp * forces(i,:) / masses(i) * dt**2
        end do 

        do i = 1, num_spheres
          velocities(i,:) = velocities(i,:) + 0.5_dp * forces(i,:) / masses(i) * dt
        end do


    end subroutine do_velocity_verlet_step

  end module velocity_verlet