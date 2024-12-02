module output_helpers
  use types_module
  implicit none

contains

  subroutine print_asterisk_row(n)
    integer, intent(in) :: n
    integer :: i
    do i = 1, n
      write(*, '(A)', advance='no') '*'
    end do
    write(*,*)  ! Move to the next line after printing the asterisks
  end subroutine print_asterisk_row


end module output_helpers
