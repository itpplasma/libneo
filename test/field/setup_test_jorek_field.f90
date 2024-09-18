program setup_test_jorek_field
use, intrinsic :: iso_fortran_env, only: dp => real64
use util_for_test, only: print_test, print_ok, print_fail
use hdf5_tools, only: hid_t, h5_init, h5_create, h5_add, h5_close, h5_deinit

implicit none


call create_homogenous_field_jorek_output()


contains


subroutine create_homogenous_field_jorek_output()
    integer, parameter :: n_var = 17, n_R = 100, n_Z = 100, n_phi = 33
    integer(hid_t) :: file_id
    character(len=*), parameter :: filename = 'test_jorek_output.h5'
    character(len=100) :: comment, description
    integer :: ndim(3)
    integer, parameter :: index_now = 40000
    real(dp), dimension(:,:,:,:), allocatable :: values
    real(dp), parameter :: t_now = 198906, time = 0.0916456
    character(len=5), parameter :: variables(17) = (/"p","h","i"," "," "," "," "," ", & 
                                                   " "," "," "," ","A","_","R"," "," "/)

    comment = 'Output produced by jorek2_postproc command "rectangular_torus"'
    description = 'Diagnostic export from JOREK'
    ndim = (/n_R, n_Z, n_phi/)
    allocate(values(n_var, n_R, n_Z, n_phi))
    
    call print_test("create_homogenous_field_jorek_output")
    call h5_init()
    call h5_create(filename, file_id)
    call h5_add(file_id, 'comment', comment)
    call h5_add(file_id, 'description', description)
    call h5_add(file_id, 'index_now', index_now)
    call h5_add(file_id, 'n_var', n_var)
    call h5_add(file_id, 'dim', ndim, lbounds=(/1/), ubounds=(/4/), &
                                      comment='quantity nR nZ nPhi', unit='')
    call h5_add(file_id, 't_now', t_now)
    call h5_add(file_id, 'time', time)
    call h5_add(file_id, 'values', values, lbounds=(/1,1,1,1/), &
                                           ubounds=(/n_var,n_R,n_Z,n_phi/), &
                                           comment='Magnetic field values', unit='T')
    call h5_add(file_id, 'variables', variables, lbounds=(/1/), ubounds=(/17/))
    call h5_close(file_id)
  
    call h5_deinit()
    
    call print_ok
end subroutine create_homogenous_field_jorek_output

end program setup_test_jorek_field
