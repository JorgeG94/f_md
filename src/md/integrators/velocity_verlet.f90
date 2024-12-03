module verlet_integrator
   use stdlib_kinds, only: dp
   use center_of_mass_utils, only: HardSphereType
   use app_output
   implicit none
   private
   public :: do_weird_verlet_step
contains

   subroutine do_weird_verlet_step(hard_spheres, forces, dt, masses, old_forces, velocities, is_first_step)
      type(HardSphereType), intent(inout) :: hard_spheres
      real(dp), intent(in) :: forces(:,:)
      real(dp), intent(in) :: dt
      real(dp), intent(in) :: masses(:)
      real(dp), intent(in) :: old_forces(:,:)
      real(dp), intent(inout) :: velocities(:,:)
      logical, intent(in) :: is_first_step
      integer :: i, num_spheres

      num_spheres = hard_spheres%num_spheres
      ! call print_array(old_forces, 'NUMPY')
      ! call print_array(forces, 'NUMPY')
      if(.not. is_first_step) then
         do i = 1, num_spheres
            velocities(i,:) = velocities(i,:) - 0.5_dp * dt * (forces(i,:) + old_forces(i,:)) / masses(i)
         end do
      endif

      do i = 1, num_spheres
         hard_spheres%coords(i,:) = hard_spheres%coords(i,:) + &
            velocities(i,:) * dt - 0.5_dp * dt * dt * forces(i,:) / masses(i)
      end do

   end subroutine do_weird_verlet_step

end module verlet_integrator
