module molecule_utils 
  use periodic_table, only: symbol_to_atomic_number, atomic_number_to_mass, atomic_number_to_symbol
  use stdlib_kinds, only: dp
  implicit none 
  private

  public :: Molecule_type

  type :: Molecule_type
    integer :: num_atoms
    character(len=2), allocatable :: atom_symbols(:)
    real(dp), allocatable :: coords(:,:)
    integer, allocatable :: atomic_numbers(:)
    real(dp), allocatable :: masses(:)
    contains
    procedure :: construct_molecule
    procedure :: print_molecule
  end type Molecule_type

  contains 

  subroutine construct_molecule(this, atom_symbols, coords)
    class(Molecule_type), intent(inout) :: this
    character(len=2), intent(in) :: atom_symbols(:)
    real(dp), intent(in) :: coords(:,:)
    integer :: i

    this%num_atoms = size(atom_symbols)
    allocate(this%atom_symbols(this%num_atoms))
    allocate(this%coords(this%num_atoms,3))
    allocate(this%atomic_numbers(this%num_atoms))
    allocate(this%masses(this%num_atoms))

    this%atom_symbols = atom_symbols
    this%coords = coords

    do i = 1, this%num_atoms
      this%atomic_numbers(i) = symbol_to_atomic_number(this%atom_symbols(i))
      this%masses(i) = atomic_number_to_mass(this%atomic_numbers(i))
    end do
  end subroutine construct_molecule

  subroutine print_molecule(molecule)
    class(Molecule_type), intent(in) :: molecule
    integer :: i

    print *, "Number of atoms:", molecule%num_atoms
    print *, "Atomic symbols:"
    do i = 1, molecule%num_atoms
      print *, molecule%atom_symbols(i)
    end do
    print *, "Atomic numbers:"
    do i = 1, molecule%num_atoms
      print *, molecule%atomic_numbers(i)
    end do
    print *, "Masses:"
    do i = 1, molecule%num_atoms
      print *, molecule%masses(i)
    end do
    print *, "Coordinates:"
    do i = 1, molecule%num_atoms
      print *, molecule%coords(i,:)
    end do
  end subroutine print_molecule

end module molecule_utils