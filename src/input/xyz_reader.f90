module xyz_reader
    implicit none
    private
    public :: read_xyz

contains

    ! Subroutine to read an XYZ file
    subroutine read_xyz(filename, num_atoms, atom_symbols, coords)
        use stdlib_kinds, only: dp
        implicit none
        ! Input
        character(len=*), intent(in) :: filename  ! Input file name

        ! Output
        integer, intent(out) :: num_atoms                  ! Number of atoms
        character(len=2), allocatable, intent(out) :: atom_symbols(:)  ! Atomic symbols
        real(dp), allocatable, intent(out) :: coords(:,:)  ! Cartesian coordinates (x, y, z)

        ! Local variables
        integer :: i, unit, ios
        character(len=256) :: line
        real(dp) :: x, y, z
        character(len=2) :: symbol

        ! Open the file
        open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print *, "Error opening file:", trim(filename)
            stop
        end if

        ! Read the number of atoms
        read(unit, '(A)', iostat=ios) line
        if (ios /= 0) then
            print *, "Error reading number of atoms"
            stop
        end if
        read(line, *) num_atoms

        ! Skip the second line (title or blank)
        read(unit, '(A)', iostat=ios)

        ! Allocate arrays
        allocate(atom_symbols(num_atoms))
        allocate(coords(num_atoms,3))

        ! Read the atomic symbols and coordinates
        do i = 1, num_atoms
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) then
                print *, "Error reading atom data at line", i+2
                stop
            end if
            read(line, *) symbol, x, y, z
            atom_symbols(i) = symbol
            coords(i,:) = [x, y, z]
        end do

        ! Close the file
        close(unit)

    end subroutine read_xyz

end module xyz_reader

