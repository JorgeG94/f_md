program main
   use xyz_reader, only: read_xyz
   use stdlib_kinds, only: dp
   use stdlib_codata
   use molecule_utils, only: Molecule_type
   use center_of_mass_utils, only: HardSphereType, FragmentType, construct_hard_spheres, print_hard_spheres
   use lennard_jones, only: calculate_lj_potential, calculate_lj_forces
   use verlet_integrator, only: do_weird_verlet_step
   use app_output
   implicit none

   character(len=256) :: filename
   integer :: num_atoms
   character(len=2), allocatable :: atom_symbols(:)
   real(dp), allocatable :: coords(:,:), old_forces(:,:)
   real(dp), allocatable :: forces(:,:), velocities(:,:)
   real(dp) :: lj_potential, epsilon, sigma, dt, num_steps, scaled_epsilon, scaled_sigma
   integer :: i, j
   real(dp) :: scale_factor, kinetic_energy, temperature
   type(Molecule_type) :: molecule
   type(FragmentType) :: water
   type(HardSphereType) :: hard_spheres
   real(dp), allocatable :: potential_energies(:), kinetic_energies(:),&
      total_energies(:), temperatures(:)
   real(dp) :: total_momentum(3)
   logical :: is_first_step = .true.
   logical :: energy_conserved = .true.
   real(dp) :: Ediff, total_mass

!   epsilon = 0.155_dp
!   sigma = 3.166_dp
   !epsilon = 0.997_dp 
   ! epsilon in J 
   epsilon = 1.5d-21
   ! sigma = nm 
   sigma = 0.346231_dp
   scaled_epsilon = epsilon / epsilon
   scaled_sigma = sigma / sigma

   dt = 0.001_dp
   num_steps = 5
   allocate(potential_energies(num_steps))
   allocate(kinetic_energies(num_steps))
   allocate(total_energies(num_steps))
   allocate(temperatures(num_steps))

   filename = 'inputs/example.xyz'
   call read_xyz(filename, num_atoms, atom_symbols, coords)
   print *, "the boltzman constant ", BOLTZMANN_CONSTANT%value

   call molecule%construct_molecule(atom_symbols, coords)



   water%num_atoms = 3
   allocate(water%atom_symbols(water%num_atoms))
   water%atom_symbols = ['O ', 'H ', 'H ']
   ! this could be made object oriented
   call construct_hard_spheres(molecule, water, hard_spheres)
   ! transform coords to nm 
   hard_spheres%coords = hard_spheres%coords * 0.1_dp
   ! transform to reduced units nm / nm  
   hard_spheres%coords = hard_spheres%coords / sigma 
   ! make the mass dimensionless too 
   hard_spheres%masses = hard_spheres%masses / hard_spheres%masses(1)

   call print_array(hard_spheres%coords, 'NUMPY')
   total_mass = sum(hard_spheres%masses)

   allocate(forces(hard_spheres%num_spheres,3))
   allocate(velocities(hard_spheres%num_spheres,3))
   forces = 0.0_dp
   call random_number(scale_factor)
   velocities = 0.0_dp ! 0.01_dp * scale_factor


   ! Initial force calculation

   allocate(old_forces(hard_spheres%num_spheres,3))

   print *, "Step: ",  " Total energy: ", " Potential energy: ",&
      " Kinetic energy: "
   do i = 1, num_steps
      ! cache old forces
      old_forces = forces
      call calculate_lj_forces(hard_spheres, scaled_epsilon, scaled_sigma, lj_potential, forces)

      call do_weird_verlet_step(hard_spheres, forces, dt, hard_spheres%masses, old_forces, velocities, is_first_step)
      do j = 1, hard_spheres%num_spheres
     !    total_momentum(:) = total_momentum(:) + hard_spheres%masses(j) * velocities(j, :)
      end do
      do j = 1, hard_spheres%num_spheres
    !velocities(j,:) = velocities(j,:) - total_momentum(:) / total_mass
      end do
      ! do j = 1, hard_spheres%num_spheres
      ! print *, "O", hard_spheres%coords(j, :)
      ! end do
      is_first_step = .false.
      ! 4. Compute kinetic energy
      total_momentum = 0.0_dp
      kinetic_energy = 0.0_dp
      do j = 1, hard_spheres%num_spheres
         kinetic_energy = kinetic_energy + 0.5_dp * hard_spheres%masses(j) &
            * sum(velocities(j, :)**2)
      end do
      potential_energies(i) = lj_potential
      kinetic_energies(i) = kinetic_energy
      temperatures(i) = (2.0_dp / 3.0_dp) * kinetic_energy / (hard_spheres%num_spheres ) 
      total_energies(i) = lj_potential + kinetic_energy
      !print *, "total momentum: ", total_momentum

      ! 5. Print diagnostics
      call print_message(to_string(i) // ", " // to_string(total_energies(i)) // ", " // &
      to_string(potential_energies(i)) // ", " // to_string(kinetic_energies(i)) &
      // ", " // to_string(temperatures(i)))
   end do

   do i = 1, num_steps - 1
      Ediff = total_energies(i+1) - total_energies(i)
      !print *, "Energy difference: ", Ediff
      if(abs(Ediff).gt.1.0d-8) then
         print *, "Energy is not conserved in step ", i, " with difference ", Ediff
         energy_conserved = .false.
         stop
      end if
   end do

   if(.not. energy_conserved) then
      print *, "Energy is not conserved"
   else
      print *, "Energy is conserved"
   end if



end program main
