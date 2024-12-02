module types_module
  use stdlib_kinds, only: sp, dp
  implicit none 

  integer, parameter :: int64 = selected_int_kind(18)
  integer, parameter :: int32 = selected_int_kind(9)
  integer, parameter :: l_qp = selected_real_kind(33)

  integer, parameter :: default_int = int32
  integer, parameter :: default_real = dp

  type, public :: DPType
    real(kind=default_real) :: value
  end type DPType

  type, public :: Int4Type 
    integer(kind=default_int) :: value
  end type Int4Type

  type, public :: Int8Type
    integer(kind=int64) :: value
  end type Int8Type

end module types_module