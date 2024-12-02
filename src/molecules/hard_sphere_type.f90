module center_of_mass_utils 
  use stdlib_kinds, only: dp
  use molecule_utils, only: Molecule_type
  implicit none
  private
  public :: FragmentType, HardSphereType, compute_center_of_mass, construct_hard_spheres, print_hard_spheres 

  type :: FragmentType
    integer :: num_atoms 
    character(len=2), allocatable :: atom_symbols(:)
  end type FragmentType

  type :: HardSphereType
    integer :: num_spheres
    real(dp), allocatable :: masses(:)
    real(dp), allocatable :: coords(:,:)
  end type HardSphereType

  contains

  function compute_center_of_mass(masses, coords) result(center_of_mass)
    real(dp), intent(in) :: masses(:)
    real(dp), intent(in) :: coords(:,:)
    real(dp) :: center_of_mass(3)
    real(dp) :: total_mass
    integer :: i

    center_of_mass = 0.0_dp
    total_mass = sum(masses)
    do i = 1, size(masses)
      center_of_mass = center_of_mass + masses(i) * coords(i,:)
    end do 
    center_of_mass = center_of_mass / total_mass
  end function compute_center_of_mass

  subroutine construct_hard_spheres(molecule, fragment, hard_spheres)
    class(Molecule_type), intent(in) :: molecule
    type(FragmentType), intent(in) :: fragment
    type(HardSphereType), intent(out) :: hard_spheres
    real(dp), allocatable :: frag_masses(:)
    real(dp), allocatable :: frag_coords(:,:)
    real(dp) :: frag_center_of_mass(3)
    integer :: i, frag_idx, atom_idx

    hard_spheres%num_spheres = molecule%num_atoms / fragment%num_atoms
    allocate(hard_spheres%masses(hard_spheres%num_spheres))
    allocate(hard_spheres%coords(hard_spheres%num_spheres,3))

    do frag_idx = 1, hard_spheres%num_spheres
      atom_idx = (frag_idx - 1) * fragment%num_atoms + 1
      allocate(frag_masses(fragment%num_atoms))
      allocate(frag_coords(fragment%num_atoms,3))

      do i = 1, fragment%num_atoms
        frag_masses(i) = molecule%masses(atom_idx)
        frag_coords(i,:) = molecule%coords(atom_idx,:)
        atom_idx = atom_idx + 1
      end do

      frag_center_of_mass = compute_center_of_mass(frag_masses, frag_coords)

      hard_spheres%masses(frag_idx) = sum(frag_masses)
      hard_spheres%coords(frag_idx,:) = frag_center_of_mass

      deallocate(frag_masses, frag_coords)

    end do

  end subroutine construct_hard_spheres


  subroutine print_hard_spheres(hard_spheres)
    type(HardSphereType), intent(in) :: hard_spheres
    integer :: i

    print *, "Number of hard spheres:", hard_spheres%num_spheres
    print *, "Masses:"
    do i = 1, hard_spheres%num_spheres
      print *, hard_spheres%masses(i)
    end do
    print *, "Coordinates:"
    do i = 1, hard_spheres%num_spheres
      print *, hard_spheres%coords(i,:)
    end do
  end subroutine print_hard_spheres
  

end module center_of_mass_utils