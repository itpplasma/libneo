program setup_test_jorek_field
use, intrinsic :: iso_fortran_env, only: dp => real64
use util_for_test, only: print_test, print_ok, print_fail

implicit none


call create_mockup_jorek_output()


contains


subroutine create_mockup_jorek_output()
    use util, only: linspace
    use util_for_test_jorek_field, only: Rmin, Rmax, Zmin, Zmax, phimin, phimax, &
                                     n_var, n_R, n_Z, n_phi, ndim, index_now, &
                                     t_now, time, variables, save_filename, &
                                     comment, description, filename_len
    use hdf5_tools, only: hid_t, h5_init, h5_create, h5_add, h5_close, h5_deinit

    character(len=filename_len) :: filename
    integer(hid_t) :: file_id
    integer :: dims(3)
    real(dp), dimension(:), allocatable :: R
    real(dp), dimension(:,:,:,:), allocatable :: values

    call print_test("create_mockup_jorek_output")

    dims = (/n_phi, n_Z, n_R/)
    allocate(R(n_R))
    call linspace(Rmin, Rmax, n_R, R)
    allocate(values(n_phi, n_Z, n_R, n_var))
    call get_homogenous_field(R, values)
    
    filename = make_filename(Rmin, Rmax, Zmin, Zmax, phimin, phimax)
    call save_filename(filename)

    call h5_init()
    call h5_create(filename, file_id)
    call h5_add(file_id, 'comment', comment)
    call h5_add(file_id, 'description', description)
    call h5_add(file_id, 'dim', dims, lbounds=(/1/), ubounds=(/3/), &
                                            comment='quantity nR nZ nPhi', unit='')
    call h5_add(file_id, 'index_now', index_now)
    call h5_add(file_id, 'n_var', n_var)
    call h5_add(file_id, 'ndim', ndim)
    call h5_add(file_id, 't_now', t_now)
    call h5_add(file_id, 'time', time)
    call h5_add(file_id, 'values', values, lbounds=(/1,1,1,1/), &
                                           ubounds=(/n_phi, n_Z, n_R, n_var/), &
                                           comment='Magnetic field values', unit='T')
    call h5_add(file_id, 'variables', variables, lbounds=(/1/), ubounds=(/17/))
    call h5_close(file_id)
  
    call h5_deinit()
    
    call print_ok
end subroutine create_mockup_jorek_output

function make_filename(Rmin, Rmax, Zmin, Zmax, phimin, phimax) result(filename)
    use iso_c_binding
    use util_for_test_jorek_field, only: filename_len
    real(dp), intent(in) :: Rmin, Rmax, Zmin, Zmax, phimin, phimax
    character(len=filename_len) :: filename

    character(len=filename_len) :: tmp_filename
    character(len=filename_len) :: current_directory
    integer :: ierr

    write(tmp_filename, '(A,F5.3,A,F5.3,A,F6.3,A,F5.3,A,F5.3,A,F5.3,A)') &
        'exprs_Rmin', Rmin, '_Rmax', Rmax, &
        '_Zmin', Zmin, '_Zmax', Zmax, &
        '_phimin', phimin, '_phimax', phimax, &
        '.h5'
    call getcwd(current_directory, ierr)
    filename = trim(adjustl(current_directory)) // '/' // trim(tmp_filename)
end function make_filename

subroutine get_homogenous_field(R, values)
    real(dp), dimension(:), intent(in) :: R
    real(dp), dimension(:,:,:,:), intent(out) :: values

    integer :: dims(4)
    real(dp), dimension(:,:,:), allocatable :: A_Z, B_phi

    values = 0.0_dp
    dims = shape(values)
    allocate(A_Z(dims(1), dims(2), dims(3)))
    allocate(B_phi(dims(1), dims(2), dims(3)))
    A_Z = 1.0_dp * spread(spread(R, dim=1, ncopies=dims(1)), dim=2, ncopies=dims(2))
    B_phi = -1.0_dp
    values(:, :, :, 2) = A_Z
    values(:, :, :, 14) = B_phi
end subroutine get_homogenous_field


end program setup_test_jorek_field
