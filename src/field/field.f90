module neo_field
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neo_field_base, only: field_t
    use neo_example_field, only: example_field_t
    use neo_biotsavart_field, only: biotsavart_field_t
    use neo_polylag_field, only: polylag_field_t
    implicit none

    contains

    subroutine create_field(field, field_type, ampl, ampl2, &
                                               coils_file, &
                                               limits, field_to_interpolate, n_points)
        class(field_t), allocatable, intent(inout) :: field
        character(*), intent(in) :: field_type
        real(dp), intent(in), optional :: ampl, ampl2
        character(*), intent(in), optional :: coils_file
        real(dp), dimension(3,2), intent(in), optional :: limits
        class(field_t), intent(in), optional :: field_to_interpolate
        integer, dimension(3), intent(in), optional :: n_points

        select case(field_type)
            case("example")
                allocate(field, source=create_example_field(ampl, ampl2))
            case("biotsavart")
                allocate(field, source=create_biotsavart_field(coils_file))
            case("polylag")
                allocate(field, source= & 
                            create_polylag_field(limits, field_to_interpolate, n_points))
            case default
                print *, "Invalid field type"
                error stop
        end select
    end subroutine create_field


    function create_example_field(ampl, ampl2) result(example_field)

        real(dp), intent(in), optional :: ampl, ampl2
        class(example_field_t), allocatable :: example_field

        allocate(example_field)
        call example_field%example_field_init(ampl, ampl2)
    end function create_example_field


    function create_biotsavart_field(coils_file) result(biotsavart_field)

        character(*), intent(in), optional :: coils_file

        class(biotsavart_field_t), allocatable :: biotsavart_field

        allocate(biotsavart_field)
        call biotsavart_field%biotsavart_field_init(coils_file)
    end function create_biotsavart_field

    function create_polylag_field(limits, field, n_points) result(polylag_field)

        real(dp), dimension(3,2), intent(in) :: limits
        class(field_t), intent(in), optional :: field
        integer, dimension(3), intent(in), optional :: n_points
        class(polylag_field_t), allocatable :: polylag_field

        allocate(polylag_field)
        call polylag_field%polylag_field_init(limits, field, n_points)
    end function create_polylag_field

end module neo_field
