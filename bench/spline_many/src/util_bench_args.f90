module util_bench_args
    use, intrinsic :: iso_fortran_env, only: int64
    implicit none
    private

    public :: get_env_int

contains

    integer function get_env_int(name, default_value) result(v)
        character(len=*), intent(in) :: name
        integer, intent(in) :: default_value

        character(len=64) :: buf
        integer :: stat, n, ios

        v = default_value
        call get_environment_variable(name, buf, length=n, status=stat)
        if (stat /= 0) return
        if (n <= 0) return

        read (buf(1:n), *, iostat=ios) v
        if (ios /= 0) then
            v = default_value
        end if
    end function get_env_int

end module util_bench_args

