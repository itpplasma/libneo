module util_for_test_jorek_field
implicit none

integer, parameter :: dp = kind(1.0d0)

real(dp), parameter :: Rmin = 1.0_dp, Rmax = 2.0_dp
real(dp), parameter :: Zmin = -1.0_dp, Zmax = 1.0_dp
real(dp), parameter :: phimin = 0.0_dp, phimax = 6.283
integer, parameter :: n_var = 17, n_R = 30, n_phi = 10, n_Z = 20
integer, parameter :: ndim = 3, index_now = 1000
real(dp), parameter :: t_now = 10000, time = 0.01
character(len=5), parameter :: variables(17) = (/"p","h","i"," "," "," "," "," ", &
                                                " "," "," "," ","A","_","R"," "," "/)
character(len=*), parameter :: save_location = '../jorek_filename.txt'
character(len=*), parameter :: comment = 'Mockup output from JOREK'
character(len=*), parameter :: description = 'Diagnostic export from JOREK'
integer, parameter :: filename_len = 512


contains

subroutine save_filename(filename)
    character(len=*), intent(in) :: filename

    integer :: unit_id

    open(newunit=unit_id, file=save_location, status='replace')
    write(unit_id, '(A)') filename
    close(unit_id)
end subroutine save_filename


subroutine get_filename(filename)
    character(len=*), intent(out) :: filename

    integer :: unit_id

    open(newunit=unit_id, file=save_location, status='old')
    read(unit_id, '(A)') filename
    close(unit_id)
end subroutine get_filename


subroutine remove_saved_filename()
    logical :: exists
    integer :: unit_id, stat

    inquire(file=save_location, exist=exists)
    if (exists) open(newunit=unit_id, iostat=stat, file=save_location, status='old')
    if (stat == 0) close(unit_id, status='delete')
end subroutine remove_saved_filename

end module util_for_test_jorek_field
