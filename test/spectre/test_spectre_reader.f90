program test_spectre_reader
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use hdf5_tools, only: HID_T, h5_init, h5_deinit, h5_create, h5_close, h5_add, &
                          h5_define_group, h5_close_group
    use spectre_reader, only: spectre_data_t, load_spectre, free_spectre
    use util_for_test, only: print_test, print_ok, print_fail

    implicit none

    real(dp), parameter :: tol = 1.0e-15_dp
    integer, parameter :: lrad_ref(3) = [8, 4, 4]
    integer, parameter :: im_ref(23) = [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, &
                                        3, 3, 3, 3, 3, 4, 4, 4, 4, 4]
    integer, parameter :: in_ref(23) = [0, 5, 10, -10, -5, 0, 5, 10, -10, -5, &
                                        0, 5, 10, -10, -5, 0, 5, 10, -10, -5, &
                                        0, 5, 10]

    character(len=1024) :: fixture
    integer :: nfail

    nfail = 0
    call get_command_argument(1, fixture)
    if (len_trim(fixture) == 0) then
        print *, 'usage: test_spectre_reader.x <spectre_test.h5>'
        error stop
    end if

    call test_load_fixture
    call test_reload_after_free
    call test_reject_non_toroidal

    if (nfail > 0) error stop

contains

    subroutine test_load_fixture
        type(spectre_data_t) :: data
        integer :: ierr, nfail_before

        call print_test('load_spectre reads sizes, modes, and coefficients')
        nfail_before = nfail

        call load_spectre(trim(fixture), data, ierr)
        call check_int('ierr', ierr, 0)
        call check_sizes(data)
        call check_profiles(data)
        call check_interface_geometry(data)
        call check_vector_potential(data)

        call free_spectre(data)
        if (allocated(data%vol)) call fail('free_spectre left vol allocated')
        if (allocated(data%Rbc)) call fail('free_spectre left Rbc allocated')

        call report(nfail_before)
    end subroutine test_load_fixture

    subroutine test_reload_after_free
        type(spectre_data_t) :: data
        integer :: ierr, nfail_before

        call print_test('second load/free cycle succeeds')
        nfail_before = nfail

        call load_spectre(trim(fixture), data, ierr)
        call check_int('first ierr', ierr, 0)
        call free_spectre(data)
        call load_spectre(trim(fixture), data, ierr)
        call check_int('second ierr', ierr, 0)
        call check_dp('v3 Ate(4,11)', data%vol(3)%Ate(4, 11), &
                      9.670688715652076e-8_dp)
        call free_spectre(data)

        call report(nfail_before)
    end subroutine test_reload_after_free

    subroutine test_reject_non_toroidal
        character(len=*), parameter :: invalid = 'spectre_invalid_geometry.h5'
        integer, parameter :: igeometry_cylindrical = 2

        type(spectre_data_t) :: data
        integer :: ierr, nfail_before, iunit
        integer(HID_T) :: h5id, grpid

        call print_test('load_spectre rejects Igeometry /= 3')
        nfail_before = nfail

        call h5_init()
        call h5_create(invalid, h5id)
        call h5_define_group(h5id, 'input', grpid)
        call h5_close_group(grpid)
        call h5_define_group(h5id, 'input/physics', grpid)
        call h5_close_group(grpid)
        call h5_add(h5id, 'input/physics/Igeometry', igeometry_cylindrical)
        call h5_close(h5id)
        call h5_deinit()

        call load_spectre(invalid, data, ierr)
        if (ierr == 0) call fail('expected nonzero ierr for Igeometry = 2')
        call check_int('Igeometry', data%Igeometry, igeometry_cylindrical)
        if (allocated(data%Lrad)) call fail('partial allocation: Lrad')
        if (allocated(data%Rbc)) call fail('partial allocation: Rbc')
        if (allocated(data%vol)) call fail('partial allocation: vol')

        open (newunit=iunit, file=invalid)
        close (iunit, status='delete')

        call report(nfail_before)
    end subroutine test_reject_non_toroidal

    subroutine check_sizes(data)
        type(spectre_data_t), intent(in) :: data

        call check_int('Igeometry', data%Igeometry, 3)
        call check_int('Nfp', data%Nfp, 5)
        call check_int('Nvol', data%Nvol, 3)
        call check_int('Mvol', data%Mvol, 3)
        call check_int('Mpol', data%Mpol, 4)
        call check_int('Ntor', data%Ntor, 2)
        call check_int('mn', data%mn, 23)
    end subroutine check_sizes

    subroutine check_profiles(data)
        type(spectre_data_t), intent(in) :: data

        if (any(data%Lrad /= lrad_ref)) call fail('Lrad mismatch')
        if (any(data%im /= im_ref)) call fail('im mismatch')
        if (any(data%in /= in_ref)) call fail('in mismatch')
        call check_dp('tflux(1)', data%tflux(1), 0.1111111111111111_dp)
        call check_dp('tflux(3)', data%tflux(3), 1.0_dp)
        call check_dp('pflux(2)', data%pflux(2), 0.09585332780449957_dp)
        call check_dp('pflux(3)', data%pflux(3), 0.2569660201543683_dp)
        call check_dp('mu(1)', data%mu(1), -2.8507863819722876e-4_dp)
        call check_dp('mu(3)', data%mu(3), 5.295082606351226e-6_dp)
        call check_dp('pressure(2)', data%pressure(2), 0.0_dp)
    end subroutine check_profiles

    subroutine check_interface_geometry(data)
        type(spectre_data_t), intent(in) :: data

        call check_int('lbound(Rbc,2)', lbound(data%Rbc, 2), 0)
        call check_int('ubound(Rbc,2)', ubound(data%Rbc, 2), 3)
        call check_int('size(Rbc,1)', size(data%Rbc, 1), 23)
        call check_dp('Rbc(1,0)', data%Rbc(1, 0), 10.000000000000004_dp)
        call check_dp('Rbc(1,3)', data%Rbc(1, 3), 10.0_dp)
        call check_dp('Rbc(6,2)', data%Rbc(6, 2), 0.6666666666666666_dp)
        call check_dp('Zbs(6,1)', data%Zbs(6, 1), -0.3333333333333333_dp)
        call check_dp('Zbs(7,2)', data%Zbs(7, 2), 0.16666666666666666_dp)
        call check_dp('Zbs(6,3)', data%Zbs(6, 3), -1.0_dp)
        call check_dp('max|Rbs|', maxval(abs(data%Rbs)), 0.0_dp)
        call check_int('ubound(Zbc,2)', ubound(data%Zbc, 2), 3)
    end subroutine check_interface_geometry

    subroutine check_vector_potential(data)
        type(spectre_data_t), intent(in) :: data

        integer :: lvol

        call check_int('size(vol)', size(data%vol), 3)
        do lvol = 1, 3
            call check_int('lbound(Ate,1)', lbound(data%vol(lvol)%Ate, 1), 0)
            call check_int('ubound(Ate,1)', ubound(data%vol(lvol)%Ate, 1), &
                           lrad_ref(lvol))
            call check_int('size(Ate,2)', size(data%vol(lvol)%Ate, 2), 23)
            call check_int('ubound(Azo,1)', ubound(data%vol(lvol)%Azo, 1), &
                           lrad_ref(lvol))
        end do
        call check_dp('v1 Aze(0,1)', data%vol(1)%Aze(0, 1), &
                      0.005167453618532632_dp)
        call check_dp('v1 Ate(2,1)', data%vol(1)%Ate(2, 1), &
                      0.02652582979161685_dp)
        call check_dp('v1 Aze(2,1)', data%vol(1)%Aze(2, 1), &
                      -0.007751175421167288_dp)
        call check_dp('v1 Ate(4,1)', data%vol(1)%Ate(4, 1), &
                      -7.706473795331895e-6_dp)
        call check_dp('v1 Aze(4,1)', data%vol(1)%Aze(4, 1), &
                      -2.5132957296278044e-5_dp)
        call check_dp('v1 Ate(8,3)', data%vol(1)%Ate(8, 3), &
                      -1.204224102290733e-12_dp)
        call check_dp('v2 Ate(1,1)', data%vol(2)%Ate(1, 1), &
                      0.053055225440943626_dp)
        call check_dp('v2 Aze(1,1)', data%vol(2)%Aze(1, 1), &
                      -0.015472225828772905_dp)
        call check_dp('v2 Ate(3,10)', data%vol(2)%Ate(3, 10), &
                      5.507019253210074e-9_dp)
        call check_dp('v3 Ate(1,1)', data%vol(3)%Ate(1, 1), &
                      0.08842575129402366_dp)
        call check_dp('v3 Aze(1,1)', data%vol(3)%Aze(1, 1), &
                      -0.025768354330511768_dp)
        call check_dp('v3 Ate(4,11)', data%vol(3)%Ate(4, 11), &
                      9.670688715652076e-8_dp)
        call check_dp('v3 Aze(4,11)', data%vol(3)%Aze(4, 11), &
                      -3.46330218588557e-7_dp)
        call check_dp('max|Ato| v1', maxval(abs(data%vol(1)%Ato)), 0.0_dp)
        call check_dp('max|Azo| v2', maxval(abs(data%vol(2)%Azo)), 0.0_dp)
        call check_dp('max|Ato| v3', maxval(abs(data%vol(3)%Ato)), 0.0_dp)
    end subroutine check_vector_potential

    subroutine check_int(name, got, want)
        character(len=*), intent(in) :: name
        integer, intent(in) :: got, want

        if (got /= want) then
            print *, '    ', name, ': got ', got, ' want ', want
            nfail = nfail + 1
        end if
    end subroutine check_int

    subroutine check_dp(name, got, want)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: got, want

        if (abs(got - want) > tol) then
            print *, '    ', name, ': got ', got, ' want ', want
            nfail = nfail + 1
        end if
    end subroutine check_dp

    subroutine fail(message)
        character(len=*), intent(in) :: message

        print *, '    ', message
        nfail = nfail + 1
    end subroutine fail

    subroutine report(nfail_before)
        integer, intent(in) :: nfail_before

        if (nfail == nfail_before) then
            call print_ok
        else
            call print_fail
        end if
    end subroutine report

end program test_spectre_reader
