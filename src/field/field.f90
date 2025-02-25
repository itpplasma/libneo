module neo_field
    use neo_field_base, only: field_t
    use neo_example_field, only: example_field_t
    use neo_biotsavart_field, only: biotsavart_field_t
    use neo_polylag_field, only: polylag_field_t
    use neo_spline_field, only: spline_field_t
    use neo_field_mesh, only: field_mesh_t
    use neo_jorek_field, only: jorek_field_t
    implicit none
    integer, parameter :: dp = kind(1.0d0)

    contains

    subroutine create_field(field, field_type, ampl, ampl2, &
                                               coils_file, &
                                               limits, field_to_interpolate, n_points, &
                                               field_mesh, &
                                               filename)
        class(field_t), allocatable, intent(inout) :: field
        character(*), intent(in) :: field_type
        real(dp), intent(in), optional :: ampl, ampl2
        character(*), intent(in), optional :: coils_file
        real(dp), dimension(3,2), intent(in), optional :: limits
        class(field_t), intent(in), optional :: field_to_interpolate
        integer, dimension(3), intent(in), optional :: n_points
        class(field_mesh_t), intent(in), optional :: field_mesh
        character(*), intent(in), optional :: filename

        select case(field_type)
            case("example")
                allocate(field, source=create_example_field(ampl, ampl2))
            case("biotsavart")
                allocate(field, source=create_biotsavart_field(coils_file))
            case("polylag")
                allocate(field, source= &
                            create_polylag_field(limits, field_to_interpolate, n_points))
            case("spline")
                if (present(filename)) then
                    allocate(field, source=create_spline_field_from_file(filename))
                else if (present(field_mesh)) then
                    allocate(field, source=create_spline_field_from_mesh(field_mesh))
                else
                    allocate(field, source= &
                            create_spline_field(limits, field_to_interpolate, n_points))
                end if
            case("jorek")
                allocate(field, source=create_jorek_field(filename))
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


    function create_spline_field(limits, field, n_points) result(spline_field)
        real(dp), dimension(3,2), intent(in) :: limits
        class(field_t), intent(in), optional :: field
        integer, dimension(3), intent(in), optional :: n_points
        class(spline_field_t), allocatable :: spline_field

        type(field_mesh_t) :: field_mesh

        call field_mesh%field_mesh_init_with_field(limits, field, n_points)
        allocate(spline_field, source=create_spline_field_from_mesh(field_mesh))
    end function create_spline_field


    function create_spline_field_from_file(filename) result(spline_field_from_file)
        character(*), intent(in) :: filename
        class(spline_field_t), allocatable :: spline_field_from_file

        type(field_mesh_t) :: field_mesh

        call field_mesh%field_mesh_init_with_file(filename)
        allocate(spline_field_from_file, source=create_spline_field_from_mesh(field_mesh))

    end function create_spline_field_from_file


    function create_spline_field_from_mesh(field_mesh) result(spline_field)
        type(field_mesh_t), intent(in) :: field_mesh
        class(spline_field_t), allocatable :: spline_field

        allocate(spline_field)
        call spline_field%spline_field_init(field_mesh)
    end function create_spline_field_from_mesh


    function create_jorek_field(filename) result(jorek_field)
        character(*), intent(in) :: filename
        class(jorek_field_t), allocatable :: jorek_field

        allocate(jorek_field)
        call jorek_field%jorek_field_init(filename)
    end function create_jorek_field

end module neo_field
