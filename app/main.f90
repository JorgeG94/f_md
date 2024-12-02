program main
  use xyz_reader, only: read_xyz
  use stdlib_kinds, only: dp
  use stdlib_codata
  use molecule_utils, only: Molecule_type
  use center_of_mass_utils, only: HardSphereType, FragmentType, construct_hard_spheres, print_hard_spheres
  use lennard_jones, only: calculate_lj_potential, calculate_lj_forces
  use velocity_verlet, only: do_velocity_verlet_step
  use app_output
  implicit none

  character(len=256) :: filename
  integer :: num_atoms 
  character(len=2), allocatable :: atom_symbols(:)
  real(dp), allocatable :: coords(:,:)
  real(dp), allocatable :: forces(:,:), velocities(:,:)
  real(dp) :: lj_potential, epsilon, sigma, dt, num_steps, scaled_epsilon, scaled_sigma 
  integer :: i, j
  real(dp) :: scale_factor, kinetic_energy, temperature
  type(Molecule_type) :: molecule
  type(FragmentType) :: water 
  type(HardSphereType) :: hard_spheres

  epsilon = 0.155_dp
  sigma = 3.166_dp
  scaled_epsilon = epsilon / epsilon
  scaled_sigma = sigma / sigma

  dt = 0.001_dp
  num_steps = 3

  filename = 'inputs/example.xyz'
  call read_xyz(filename, num_atoms, atom_symbols, coords)
  print *, "the boltzman constant ", BOLTZMANN_CONSTANT%value

  call molecule%construct_molecule(atom_symbols, coords)



  water%num_atoms = 3
  allocate(water%atom_symbols(water%num_atoms))
  water%atom_symbols = ['O ', 'H ', 'H ']
  ! this could be made object oriented
  call construct_hard_spheres(molecule, water, hard_spheres)
    print *, "Positions:"
    call print_array(hard_spheres%coords, 'NUMPY')
  hard_spheres%masses = hard_spheres%masses / hard_spheres%masses(1)
  hard_spheres%coords = hard_spheres%coords / sigma

  allocate(forces(hard_spheres%num_spheres,3))
  allocate(velocities(hard_spheres%num_spheres,3))
  forces = 0.0_dp
  call random_number(scale_factor)
  velocities = 0.01_dp * scale_factor

 
  ! Initial force calculation
call calculate_lj_forces(hard_spheres, epsilon, sigma, lj_potential, forces)


do i = 1, num_steps
    print *, "Step:", i

    ! 1. Velocity Verlet Step: Update positions and first half-step of velocities
    call do_velocity_verlet_step(hard_spheres, forces, dt, hard_spheres%masses, velocities)

    ! 2. Recalculate forces based on updated positions
    call calculate_lj_forces(hard_spheres, epsilon, sigma, lj_potential, forces)

    ! 3. Complete the second half-step of velocities
    do j = 1, hard_spheres%num_spheres
        velocities(j,:) = velocities(j,:) + 0.5_dp * forces(j,:) / hard_spheres%masses(j) * dt
    end do

    ! 4. Compute kinetic energy
    kinetic_energy = 0.0_dp
    do j = 1, hard_spheres%num_spheres
        kinetic_energy = kinetic_energy + 0.5_dp * hard_spheres%masses(j) * sum(velocities(j, :)**2)
    end do

    ! 5. Print diagnostics
    print *, "Positions:"
    call print_array(hard_spheres%coords, 'NUMPY')
    print *, "Potential energy: ", lj_potential
    print *, "Kinetic energy: ", kinetic_energy
    print *, "Total energy: ", lj_potential + kinetic_energy
end do








end program main
