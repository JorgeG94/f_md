module first_steps
  implicit none
  private

  public :: say_hello
contains
  subroutine say_hello
    print *, "Hello, first_steps!"
  end subroutine say_hello
end module first_steps
