module periodic_table 
  use stdlib_kinds, only: dp
  implicit none 
  private 
  public :: symbol_to_atomic_number, atomic_number_to_mass, atomic_number_to_symbol

  integer, parameter :: max_elements = 17

  character(len=2), dimension(max_elements) :: symbols = [ &
    'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl' &
  ]
  integer, dimension(max_elements) :: atomic_numbers = [ &
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17 &
  ]
  real(dp), dimension(max_elements) :: masses = [ &
    1.00794, 4.002602, 6.941, 9.012182, 10.811, 12.0107, 14.0067, 15.9994, 18.9984032,&
    20.1797, 22.98976928, 24.3050, 26.9815386, 28.0855, 30.973762, 32.065, 35.453 &
  ]

  contains 
  
  function symbol_to_atomic_number(symbol) result(atomic_number)
    character(len=*), intent(in) :: symbol
    integer :: atomic_number
    integer :: i

    atomic_number = -1
    do i = 1, max_elements
      if (trim(symbol) == trim(symbols(i))) then
        atomic_number = atomic_numbers(i)
        return 
      end if
    end do
  end function symbol_to_atomic_number

  function atomic_number_to_mass(atomic_number) result(mass)
    integer, intent(in) :: atomic_number
    real(dp) :: mass
    integer :: i

    mass = -1.0_dp
    do i = 1, max_elements
      if (atomic_number == atomic_numbers(i)) then
        mass = masses(i)
        return
      end if
    end do
  end function atomic_number_to_mass

  function atomic_number_to_symbol(atomic_number) result(symbol)
    integer, intent(in) :: atomic_number
    character(len=2) :: symbol
    integer :: i

    symbol = '??'
    do i = 1, max_elements
      if (atomic_number == atomic_numbers(i)) then
        symbol = symbols(i)
        return
      end if
    end do
  end function atomic_number_to_symbol

end module periodic_table