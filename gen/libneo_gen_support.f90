! Shared leaf-writing support for libneo FortSym generators.
! One emit path and one FortSym revision source: LIBNEO_FORTSYM_REVISION in
! CMakeLists.txt, delivered as the configured include file
! libneo_fortsym_revision.inc.
module libneo_gen_support
    use fortsym_expr, only: expr_t
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: native_engine_t, make_native_engine
    use fortsym_string, only: str_t, str, chars
    use fortsym_kernel, only: kernel_spec_t, emit_kernel
    use fortgen_horner_emit, only: horner_kernel_spec_t, emit_horner_kernel
    use fortgen_derivative_set, only: partials, radial_mixed
    use fortgen_kernel_target, only: TARGET_FORTRAN_OPENACC
    use fortgen_string, only: fg_str_t => str_t, fortgen_str => str, &
                              fortgen_chars => chars
    implicit none
    private
    include 'libneo_fortsym_revision.inc'
    public :: fortsym_revision, leaf_writer_t, write_spline_leaf

    !> Writes generated leaves into one file. The writer owns the output unit,
    !> the simplifying engine of the current FortSym context and the leaf
    !> conventions; kernel inputs bound with fortsym_jet's bind apply.
    type :: leaf_writer_t
        type(native_engine_t) :: engine
        character(:), allocatable :: generator
        integer :: unit = -1
    contains
        procedure :: open => writer_open
        procedure :: line => writer_line
        procedure :: close => writer_close
        procedure :: emit => writer_emit
        procedure :: expand => writer_expand
        procedure :: simplify => writer_simplify
    end type
contains
    subroutine writer_open(w, generator, path)
        class(leaf_writer_t), intent(inout) :: w
        character(*), intent(in) :: generator, path

        w%engine = make_native_engine()
        w%generator = generator
        open (newunit=w%unit, file=path, status='replace', action='write')
    end subroutine

    subroutine writer_line(w, text)
        class(leaf_writer_t), intent(in) :: w
        character(*), intent(in) :: text
        write (w%unit, '(a)') text
    end subroutine

    subroutine writer_close(w)
        class(leaf_writer_t), intent(inout) :: w
        close (w%unit)
        w%unit = -1
    end subroutine

    !> Fully expand an expression; fails closed when the engine refuses.
    function writer_expand(w, e) result(r)
        class(leaf_writer_t), intent(inout) :: w
        type(expr_t), intent(in) :: e
        type(expr_t) :: r
        type(engine_result_t) :: res

        res = w%engine%expand(e)
        if (.not. res%ok) error stop chars(res%message)
        r = res%value
    end function

    function writer_simplify(w, e) result(r)
        class(leaf_writer_t), intent(inout) :: w
        type(expr_t), intent(in) :: e
        type(expr_t) :: r
        type(engine_result_t) :: res

        res = w%engine%simplify(e)
        if (.not. res%ok) error stop chars(res%message)
        r = res%value
    end function

    !> Emit one pure elemental leaf subroutine without a module wrapper, so the
    !> consumer can include it next to its hot caller.
    subroutine writer_emit(w, name, args, outputs, roots)
        class(leaf_writer_t), intent(inout) :: w
        character(*), intent(in) :: name, args(:), outputs(:)
        type(expr_t), intent(in) :: roots(:)
        type(kernel_spec_t) :: spec
        type(str_t) :: source
        logical :: ok
        character(:), allocatable :: message
        integer :: j

        spec%name = str(name)
        spec%generator = str(w%generator)
        spec%generator_revision = str(fortsym_revision)
        spec%regenerate_command = str('cmake --build <build> --target '// &
                                      'regenerate_'//w%generator)
        spec%pure_procedure = .true.
        spec%elemental_procedure = .true.
        allocate (spec%args(size(args)), spec%outputs(size(outputs)))
        do j = 1, size(args)
            spec%args(j) = str(trim(args(j)))
        end do
        do j = 1, size(outputs)
            spec%outputs(j) = str(trim(outputs(j)))
        end do
        source = emit_kernel(roots, spec, ok, message=message)
        if (.not. ok) error stop message
        write (w%unit, '(a)') chars(source)
    end subroutine

    !> Emit the cell kernel spline<d>d_o<degree...>_<d<order>|rmix>_nq<nq|n> of
    !> a fixed-degree batch spline: value and derivatives up to `order` of the
    !> local polynomial for nq quantities (nq = 0: runtime count). `rmix` asks
    !> for the value, the gradient and the second derivatives along direction 1.
    subroutine write_spline_leaf(unit, d, degree, nq, order, rmix)
        integer, intent(in) :: unit, d, degree, nq, order
        logical, intent(in) :: rmix
        type(horner_kernel_spec_t) :: spec
        type(fg_str_t) :: source
        character(:), allocatable :: message, name
        character(16) :: text
        logical :: ok

        write (text, '(i0)') d
        name = 'spline'//trim(text)//'d_o'
        write (text, '(i0)') degree
        name = name//repeat(trim(text), d)//'_'
        if (rmix) then
            name = name//'rmix'
            spec%derivatives = radial_mixed(d)
        else
            write (text, '(i0)') order
            name = name//'d'//trim(text)
            spec%derivatives = partials(d, order)
        end if
        if (nq == 0) then
            name = name//'_nqn'
        else
            write (text, '(i0)') nq
            name = name//'_nq'//trim(text)
        end if
        spec%name = fortgen_str(name)
        spec%degree = 0
        spec%degree(1:d) = degree
        spec%nq = nq
        spec%target = TARGET_FORTRAN_OPENACC
        source = emit_horner_kernel(spec, ok, message)
        if (.not. ok) error stop message
        write (unit, '(a)') ''
        write (unit, '(a)') fortgen_chars(source)
    end subroutine
end module libneo_gen_support
