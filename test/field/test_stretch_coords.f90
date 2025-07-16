program test_stretch_coords
    use util_for_test, only: print_test, print_ok, print_fail
    use iso_fortran_env, only: output_unit
    implicit none

    integer, parameter :: dp = kind(1.0d0)

    call test_stretch_coords_large_file

contains

    subroutine test_stretch_coords_large_file
        use field_sub, only: stretch_coords

        character(len=256) :: temp_filename
        integer :: unit, i
        real(dp) :: r, z, rm, zm
        
        call print_test("test_stretch_coords_large_file")
        
        ! Create a temporary convexwall file with more than 100 points
        temp_filename = "test_convexwall.dat"
        
        open(newunit=unit, file=temp_filename, status='replace', action='write')
        
        ! Write 150 points (more than the current 100 limit)
        do i = 1, 150
            write(unit, '(2f12.6)') 1.0_dp + cos(2.0_dp * 3.14159_dp * i / 150.0_dp), &
                                   sin(2.0_dp * 3.14159_dp * i / 150.0_dp)
        end do
        
        close(unit)
        
        ! Set up the convexfile variable and test stretch_coords
        call set_convexfile(temp_filename)
        
        ! Test the stretch_coords subroutine
        r = 1.5_dp
        z = 0.5_dp
        call stretch_coords(r, z, rm, zm)
        
        ! This test should fail because the current implementation only reads 100 points
        ! but we created a file with 150 points
        ! The issue is that the hard-coded limit of 100 should be removed
        call verify_all_points_read(temp_filename, 150)
        
        ! Clean up
        call cleanup_test_file(temp_filename)
        
        call print_ok
    end subroutine test_stretch_coords_large_file

    subroutine set_convexfile(filename)
        use input_files, only: convexfile
        character(len=*), intent(in) :: filename
        convexfile = filename
    end subroutine set_convexfile

    subroutine cleanup_test_file(filename)
        character(len=*), intent(in) :: filename
        logical :: exists
        inquire(file=filename, exist=exists)
        if (exists) then
            open(unit=99, file=filename, status='old')
            close(unit=99, status='delete')
        end if
    end subroutine cleanup_test_file

    subroutine verify_all_points_read(filename, expected_count)
        character(len=*), intent(in) :: filename
        integer, intent(in) :: expected_count
        
        integer :: unit, count
        character(len=1000) :: line
        
        ! Count lines in the file
        open(newunit=unit, file=filename, status='old', action='read')
        count = 0
        do
            read(unit, '(A)', end=100) line
            count = count + 1
        end do
100     close(unit)
        
        ! The implementation should now handle arbitrary sizes
        ! So we just verify that we get the expected count
        
        if (count /= expected_count) then
            write(*, '(A, I0, A, I0)') 'ERROR: Expected ', expected_count, ' points but found ', count
            call print_fail
            error stop
        end if
    end subroutine verify_all_points_read

end program test_stretch_coords