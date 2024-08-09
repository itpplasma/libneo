module libneo_field
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use libneo_field_base, only: field_t
    use libneo_example_field, only: example_field_t
    use libneo_biotsavart_field, only: biotsavart_field_t


    implicit none


    contains


    function create_field(field_type, ampl, ampl2, coils_file) result(field)
        class(field_t), allocatable :: field

        character(*), intent(in) :: field_type
        real(dp), intent(in), optional :: ampl, ampl2
        character(*), intent(in), optional :: coils_file

        select case(field_type)
            case("example")
                field = create_example_field(ampl, ampl2)
            case("biotsavart")
                field = create_biotsavart_field(coils_file)
            case default
                print *, "Invalid field type"
                error stop
        end select
    end function create_field


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


end module libneo_field