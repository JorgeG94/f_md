module verlet_integrator
   use stdlib_kinds, only: dp
   use center_of_mass_utils, only: HardSphereType
   implicit none
   private
   public :: do_verlet_step
contains

   subroutine do_verlet_step(hard_spheres, forces, dt, masses, velocities, is_first_step)
      type(HardSphereType), intent(inout) :: hard_spheres
      real(dp), intent(in) :: forces(:,:)
      real(dp), intent(in) :: dt
      real(dp), intent(in) :: masses(:)
      real(dp), intent(inout) :: velocities(:,:)
      logical, intent(in) :: is_first_step
      integer :: i, num_spheres

      num_spheres = hard_spheres%num_spheres


      if(.not. is_first_step) then
         do i = 1, num_spheres
            velocities(i,:) = velocities(i,:) - 0.5_dp * dt * forces(i,:) / masses(i)
         end do
      endif

      do i = 1, num_spheres
         hard_spheres%coords(i,:) = hard_spheres%coords(i,:) + &
            velocities(i,:) * dt + 0.5_dp * forces(i,:) / masses(i) * dt**2
      end do

   end subroutine do_verlet_step

end module verlet_integrator
