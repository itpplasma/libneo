module test_util
    ! TODO move to libneo
    implicit none

    contains

    subroutine print_test(test_name)
        character(*) :: test_name
        print *, "==> ", test_name
    end subroutine print_test


    subroutine print_ok
        print *, "    .................................................... OK"
    end subroutine print_ok


    subroutine print_fail
        print *, "    .................................................... FAIL"
    end subroutine print_fail

end module test_util
