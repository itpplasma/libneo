program test_spectre_chart_validity
    !> The SPEC interface blend must not fold: inside a folded volume the
    !> Jacobian of (rho_g, theta, zeta) -> x changes sign and every field
    !> quantity divided by sqrt(g) is meaningless. validate_spectre_chart
    !> reports the first folded volume.
    !>
    !> Oracles: a circular two-volume tokamak (valid), the same boundary with
    !> the elongated inner interface of a degenerate force minimization whose
    !> blend folds near theta = pi/2 and 3 pi/2 (volume 2 folded), and the
    !> committed stellarator fixture (valid).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use spectre_reader, only: spectre_data_t, load_spectre, free_spectre
    use libneo_coordinates, only: coordinate_system_t, &
                                  spectre_coordinate_system_t, &
                                  make_spectre_coordinate_system, &
                                  validate_spectre_chart
    use field_spectre, only: spectre_field_t, create_spectre_field
    use util_for_test, only: print_test, print_ok, print_fail

    implicit none

    real(dp), parameter :: R_FOLDED(5) = [10.0986_dp, 0.259398_dp, &
        -0.00459824_dp, 0.00044414_dp, -0.000174349_dp]
    real(dp), parameter :: Z_FOLDED(5) = [0.0_dp, -0.982564_dp, 0.00586076_dp, &
        0.0119019_dp, 0.00015006_dp]
    real(dp), parameter :: R_CIRCLE(5) = [10.0_dp, 0.5_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    real(dp), parameter :: Z_CIRCLE(5) = [0.0_dp, -0.5_dp, 0.0_dp, 0.0_dp, 0.0_dp]

    character(len=1024) :: fixture
    integer :: nfail

    nfail = 0
    call get_command_argument(1, fixture)
    if (len_trim(fixture) == 0) then
        print *, 'usage: test_spectre_chart_validity.x <spectre.h5>'
        error stop
    end if

    call check_tokamak('circular inner interface is a valid chart', &
                       R_CIRCLE, Z_CIRCLE, 0)
    call check_tokamak('elongated inner interface folds volume 2', &
                       R_FOLDED, Z_FOLDED, 2)
    call check_fixture

    if (nfail > 0) error stop

contains

    subroutine check_tokamak(label, r_inner, z_inner, expected)
        character(*), intent(in) :: label
        real(dp), intent(in) :: r_inner(5), z_inner(5)
        integer, intent(in) :: expected

        type(spectre_data_t) :: data
        class(coordinate_system_t), allocatable :: cs
        integer :: ierr

        call print_test(label)
        data%Mvol = 2
        data%mn = 5
        allocate (data%im(5), data%in(5))
        data%im = [0, 1, 2, 3, 4]
        data%in = 0
        allocate (data%Rbc(5, 0:2), data%Zbs(5, 0:2), data%Rbs(5, 0:2), &
                  data%Zbc(5, 0:2))
        data%Rbs = 0.0_dp
        data%Zbc = 0.0_dp
        data%Rbc(:, 0) = [10.0962_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
        data%Zbs(:, 0) = 0.0_dp
        data%Rbc(:, 1) = r_inner
        data%Zbs(:, 1) = z_inner
        data%Rbc(:, 2) = [10.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
        data%Zbs(:, 2) = [0.0_dp, -1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]

        call make_spectre_coordinate_system(cs, data)
        select type (cs)
        type is (spectre_coordinate_system_t)
            call validate_spectre_chart(cs, ierr)
        end select
        if (ierr == expected) then
            call print_ok
        else
            call print_fail
            print *, '  validate_spectre_chart ierr =', ierr, ' expected', expected
            nfail = nfail + 1
        end if
    end subroutine check_tokamak

    subroutine check_fixture
        type(spectre_field_t) :: field
        integer :: ierr

        call print_test('stellarator fixture loads with a valid chart')
        call create_spectre_field(field, trim(fixture), ierr)
        if (ierr == 0) then
            call print_ok
        else
            call print_fail
            print *, '  create_spectre_field ierr =', ierr
            nfail = nfail + 1
        end if
    end subroutine check_fixture

end program test_spectre_chart_validity
