module openacc
    !> Stub module providing OpenACC API when compiling without OpenACC support.
    !> This allows the code to compile with standard Fortran compilers.
    !> When _OPENACC is defined, the real openacc module from the compiler is used.
    implicit none
    private

    public :: acc_is_present

contains

    logical function acc_is_present(a, len)
        !> Stub: always returns false when OpenACC is not available.
        class(*), intent(in) :: a
        integer, intent(in) :: len

        acc_is_present = .false.
    end function acc_is_present

end module openacc
