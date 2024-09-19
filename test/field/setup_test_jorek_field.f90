program setup_test_jorek_field
use, intrinsic :: iso_fortran_env, only: dp => real64
use util_for_test, only: print_test, print_ok, print_fail
use hdf5_tools, only: hid_t, h5_init, h5_create, h5_add, h5_close, h5_deinit

implicit none


call create_mockup_jorek_output()


contains


subroutine create_mockup_jorek_output()
    integer, parameter :: n_var = 17, n_R = 100, n_Z = 100, n_phi = 33
    integer(hid_t) :: file_id
    character(len=*), parameter :: filename = 'jorek_output.h5'
    character(len=100) :: comment, description
    integer :: dims(3)
    integer, parameter :: ndim = 3, index_now = 1000
    real(dp), dimension(:,:,:,:), allocatable :: values
    real(dp), parameter :: t_now = 10000, time = 0.01
    character(len=5), parameter :: variables(17) = (/"p","h","i"," "," "," "," "," ", & 
                                                   " "," "," "," ","A","_","R"," "," "/)

    comment = 'Output produced by jorek2_postproc command "rectangular_torus"'
    description = 'Diagnostic export from JOREK'
    dims = (/n_phi, n_Z, n_R/)
    allocate(values(n_phi, n_Z, n_R, n_var))
    call get_homogenous_field(values)
    
    call print_test("create_mockup_jorek_output")
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

subroutine get_homogenous_field(values)
    real(dp), dimension(:,:,:,:), intent(out) :: values

    integer :: dims(4), i
    real(dp), dimension(:), allocatable :: R
    real(dp), parameter :: R_min = 1.0_dp, R_max = 2.0_dp
    real(dp), dimension(:,:,:), allocatable :: Aphi, A3, Bz

    values = 0.0_dp
    dims = shape(values)
    allocate(Aphi(dims(1), dims(2), dims(3)))
    allocate(A3(dims(1), dims(2), dims(3)))
    allocate(Bz(dims(1), dims(2), dims(3)))
    Aphi = 0.5_dp
    allocate(R(dims(3)))
    R = [(R_min + (R_max - R_min)*(i - 1)/(dims(3) - 1), i = 1, dims(3))]
    A3 = Aphi * spread(spread(R, dim=1, ncopies=dims(1)), dim=2, ncopies=dims(2))
    Bz = 1.0_dp
    values(:, :, :, 3) = A3
    values(:, :, :, 13) = Bz
end subroutine get_homogenous_field


end program setup_test_jorek_field
