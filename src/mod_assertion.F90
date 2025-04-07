module mod_assertion

    use, intrinsic :: iso_fortran_env
    implicit none

    private

    ! Public methods
    public :: assert_true
    public :: assert_false
    public :: assert_eq
    public :: assert_neq

    interface assert_eq
        procedure assert_eq_single_int8
        procedure assert_eq_single_int16
        procedure assert_eq_single_int32
        procedure assert_eq_single_int64
        procedure assert_eq_single_real32
        procedure assert_eq_single_real64
        procedure assert_eq_single_real128
        procedure assert_eq_single_logical
        procedure assert_eq_single_character
        procedure assert_eq_array_int8
        procedure assert_eq_array_int16
        procedure assert_eq_array_int32
        procedure assert_eq_array_int64
        procedure assert_eq_array_real32
        procedure assert_eq_array_real64
        procedure assert_eq_array_real128
        procedure assert_eq_array_logical
        procedure assert_eq_array_character
    end interface assert_eq

    interface assert_neq
        procedure assert_neq_single_int8
        procedure assert_neq_single_int16
        procedure assert_neq_single_int32
        procedure assert_neq_single_int64
        procedure assert_neq_single_real32
        procedure assert_neq_single_real64
        procedure assert_neq_single_real128
        procedure assert_neq_single_logical
        procedure assert_neq_single_character
        procedure assert_neq_array_int8
        procedure assert_neq_array_int16
        procedure assert_neq_array_int32
        procedure assert_neq_array_int64
        procedure assert_neq_array_real32
        procedure assert_neq_array_real64
        procedure assert_neq_array_real128
        procedure assert_neq_array_logical
        procedure assert_neq_array_character
    end interface assert_neq

contains

    subroutine assert_true(val)

        ! Subroutine arguments
        logical, intent(in) :: val

        if (.not. val) error stop

    end subroutine assert_true

    subroutine assert_false(val)

        ! Subroutine arguments
        logical, intent(in) :: val

        if (val) error stop

    end subroutine assert_false

    subroutine assert_eq_single_int8(lhs, rhs)

        ! Subroutine arguments
        integer(int8), intent(in) :: lhs
        integer(int8), intent(in) :: rhs

        if (lhs /= rhs) error stop

    end subroutine assert_eq_single_int8

    subroutine assert_eq_single_int16(lhs, rhs)

        ! Subroutine arguments
        integer(int16), intent(in) :: lhs
        integer(int16), intent(in) :: rhs

        if (lhs /= rhs) error stop

    end subroutine assert_eq_single_int16

    subroutine assert_eq_single_int32(lhs, rhs)

        ! Subroutine arguments
        integer(int32), intent(in) :: lhs
        integer(int32), intent(in) :: rhs

        if (lhs /= rhs) error stop

    end subroutine assert_eq_single_int32

    subroutine assert_eq_single_int64(lhs, rhs)

        ! Subroutine arguments
        integer(int64), intent(in) :: lhs
        integer(int64), intent(in) :: rhs

        if (lhs /= rhs) error stop

    end subroutine assert_eq_single_int64

    subroutine assert_eq_single_real32(lhs, rhs)

        ! Subroutine arguments
        real(real32), intent(in) :: lhs
        real(real32), intent(in) :: rhs

        if (abs(lhs - rhs) > epsilon(1.0_real32)) error stop

    end subroutine assert_eq_single_real32

    subroutine assert_eq_single_real64(lhs, rhs)

        ! Subroutine arguments
        real(real64), intent(in) :: lhs
        real(real64), intent(in) :: rhs

        if (abs(lhs - rhs) > epsilon(1.0_real64)) error stop

    end subroutine assert_eq_single_real64

    subroutine assert_eq_single_real128(lhs, rhs)

        ! Subroutine arguments
        real(real128), intent(in) :: lhs
        real(real128), intent(in) :: rhs

        if (abs(lhs - rhs) > epsilon(1.0_real128)) error stop

    end subroutine assert_eq_single_real128

    subroutine assert_eq_single_logical(lhs, rhs)

        ! Subroutine arguments
        logical, intent(in) :: lhs
        logical, intent(in) :: rhs

        if (lhs .neqv. rhs) error stop

    end subroutine assert_eq_single_logical

    subroutine assert_eq_single_character(lhs, rhs)

        ! Subroutine arguments
        character(len=*), intent(in) :: lhs
        character(len=*), intent(in) :: rhs

        if (lhs /= rhs) error stop

    end subroutine assert_eq_single_character

    subroutine assert_eq_array_int8(lhs, rhs)

        ! Subroutine arguments
        integer(int8), intent(in) :: lhs(:)
        integer(int8), intent(in) :: rhs(:)

        if (any(lhs /= rhs)) error stop

    end subroutine assert_eq_array_int8

    subroutine assert_eq_array_int16(lhs, rhs)

        ! Subroutine arguments
        integer(int16), intent(in) :: lhs(:)
        integer(int16), intent(in) :: rhs(:)

        if (any(lhs /= rhs)) error stop

    end subroutine assert_eq_array_int16

    subroutine assert_eq_array_int32(lhs, rhs)

        ! Subroutine arguments
        integer(int32), intent(in) :: lhs(:)
        integer(int32), intent(in) :: rhs(:)

        if (any(lhs /= rhs)) error stop

    end subroutine assert_eq_array_int32

    subroutine assert_eq_array_int64(lhs, rhs)

        ! Subroutine arguments
        integer(int64), intent(in) :: lhs(:)
        integer(int64), intent(in) :: rhs(:)

        if (any(lhs /= rhs)) error stop

    end subroutine assert_eq_array_int64

    subroutine assert_eq_array_real32(lhs, rhs)

        ! Subroutine arguments
        real(real32), intent(in) :: lhs(:)
        real(real32), intent(in) :: rhs(:)

        if (any(abs(lhs - rhs) > epsilon(1.0_real32))) error stop

    end subroutine assert_eq_array_real32

    subroutine assert_eq_array_real64(lhs, rhs)

        ! Subroutine arguments
        real(real64), intent(in) :: lhs(:)
        real(real64), intent(in) :: rhs(:)

        if (any(abs(lhs - rhs) > epsilon(1.0_real64))) error stop

    end subroutine assert_eq_array_real64

    subroutine assert_eq_array_real128(lhs, rhs)

        ! Subroutine arguments
        real(real128), intent(in) :: lhs(:)
        real(real128), intent(in) :: rhs(:)

        if (any(abs(lhs - rhs) > epsilon(1.0_real128))) error stop

    end subroutine assert_eq_array_real128

    subroutine assert_eq_array_logical(lhs, rhs)

        ! Subroutine arguments
        logical, intent(in) :: lhs(:)
        logical, intent(in) :: rhs(:)

        if (any(lhs .neqv. rhs)) error stop

    end subroutine assert_eq_array_logical

    subroutine assert_eq_array_character(lhs, rhs)

        ! Subroutine arguments
        character(len=*), intent(in) :: lhs(:)
        character(len=*), intent(in) :: rhs(:)

        if (any(lhs /= rhs)) error stop

    end subroutine assert_eq_array_character

    subroutine assert_neq_single_int8(lhs, rhs)

        ! Subroutine arguments
        integer(int8), intent(in) :: lhs
        integer(int8), intent(in) :: rhs

        if (lhs == rhs) error stop

    end subroutine assert_neq_single_int8

    subroutine assert_neq_single_int16(lhs, rhs)

        ! Subroutine arguments
        integer(int16), intent(in) :: lhs
        integer(int16), intent(in) :: rhs

        if (lhs == rhs) error stop

    end subroutine assert_neq_single_int16

    subroutine assert_neq_single_int32(lhs, rhs)

        ! Subroutine arguments
        integer(int32), intent(in) :: lhs
        integer(int32), intent(in) :: rhs

        if (lhs == rhs) error stop

    end subroutine assert_neq_single_int32

    subroutine assert_neq_single_int64(lhs, rhs)

        ! Subroutine arguments
        integer(int64), intent(in) :: lhs
        integer(int64), intent(in) :: rhs

        if (lhs == rhs) error stop

    end subroutine assert_neq_single_int64

    subroutine assert_neq_single_real32(lhs, rhs)

        ! Subroutine arguments
        real(real32), intent(in) :: lhs
        real(real32), intent(in) :: rhs

        if (abs(lhs - rhs) <= epsilon(1.0_real32)) error stop

    end subroutine assert_neq_single_real32

    subroutine assert_neq_single_real64(lhs, rhs)

        ! Subroutine arguments
        real(real64), intent(in) :: lhs
        real(real64), intent(in) :: rhs

        if (abs(lhs - rhs) <= epsilon(1.0_real64)) error stop

    end subroutine assert_neq_single_real64

    subroutine assert_neq_single_real128(lhs, rhs)

        ! Subroutine arguments
        real(real128), intent(in) :: lhs
        real(real128), intent(in) :: rhs

        if (abs(lhs - rhs) <= epsilon(1.0_real128)) error stop

    end subroutine assert_neq_single_real128

    subroutine assert_neq_single_logical(lhs, rhs)

        ! Subroutine arguments
        logical, intent(in) :: lhs
        logical, intent(in) :: rhs

        if (lhs .eqv. rhs) error stop

    end subroutine assert_neq_single_logical

    subroutine assert_neq_single_character(lhs, rhs)

        ! Subroutine arguments
        character(len=*), intent(in) :: lhs
        character(len=*), intent(in) :: rhs

        if (lhs == rhs) error stop

    end subroutine assert_neq_single_character

    subroutine assert_neq_array_int8(lhs, rhs)

        ! Subroutine arguments
        integer(int8), intent(in) :: lhs(:)
        integer(int8), intent(in) :: rhs(:)

        if (all(lhs == rhs)) error stop

    end subroutine assert_neq_array_int8

    subroutine assert_neq_array_int16(lhs, rhs)

        ! Subroutine arguments
        integer(int16), intent(in) :: lhs(:)
        integer(int16), intent(in) :: rhs(:)

        if (all(lhs == rhs)) error stop

    end subroutine assert_neq_array_int16

    subroutine assert_neq_array_int32(lhs, rhs)

        ! Subroutine arguments
        integer(int32), intent(in) :: lhs(:)
        integer(int32), intent(in) :: rhs(:)

        if (all(lhs == rhs)) error stop

    end subroutine assert_neq_array_int32

    subroutine assert_neq_array_int64(lhs, rhs)

        ! Subroutine arguments
        integer(int64), intent(in) :: lhs(:)
        integer(int64), intent(in) :: rhs(:)

        if (all(lhs == rhs)) error stop

    end subroutine assert_neq_array_int64

    subroutine assert_neq_array_real32(lhs, rhs)

        ! Subroutine arguments
        real(real32), intent(in) :: lhs(:)
        real(real32), intent(in) :: rhs(:)

        if (all(abs(lhs - rhs) <= epsilon(1.0_real32))) error stop

    end subroutine assert_neq_array_real32

    subroutine assert_neq_array_real64(lhs, rhs)

        ! Subroutine arguments
        real(real64), intent(in) :: lhs(:)
        real(real64), intent(in) :: rhs(:)

        if (all(abs(lhs - rhs) <= epsilon(1.0_real64))) error stop

    end subroutine assert_neq_array_real64

    subroutine assert_neq_array_real128(lhs, rhs)

        ! Subroutine arguments
        real(real128), intent(in) :: lhs(:)
        real(real128), intent(in) :: rhs(:)

        if (all(abs(lhs - rhs) <= epsilon(1.0_real128))) error stop

    end subroutine assert_neq_array_real128

    subroutine assert_neq_array_logical(lhs, rhs)

        ! Subroutine arguments
        logical, intent(in) :: lhs(:)
        logical, intent(in) :: rhs(:)

        if (all(lhs .eqv. rhs)) error stop

    end subroutine assert_neq_array_logical

    subroutine assert_neq_array_character(lhs, rhs)

        ! Subroutine arguments
        character(len=*), intent(in) :: lhs(:)
        character(len=*), intent(in) :: rhs(:)

        if (all(lhs == rhs)) error stop

    end subroutine assert_neq_array_character

end module mod_assertion
