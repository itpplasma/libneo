! Shared leaf-writing support for libneo FortSym generators.
! One emit path and one FortSym revision source: LIBNEO_FORTSYM_REVISION in
! CMakeLists.txt, delivered as the configured include file
! libneo_fortsym_revision.inc.
module libneo_gen_support
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: native_engine_t, make_native_engine
    use fortsym_string, only: str_t, str, chars
    use fortsym_kernel, only: kernel_spec_t, kernel_binding_t, emit_kernel
    implicit none
    private
    include 'libneo_fortsym_revision.inc'
    public :: fortsym_revision, leaf_writer_t, binding

    !> Writes generated leaves into one file. The generator owns the arena; the
    !> writer owns the output unit, the simplifying engine and leaf conventions.
    type :: leaf_writer_t
        type(arena_t), pointer :: arena => null()
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
    subroutine writer_open(w, arena, generator, path)
        class(leaf_writer_t), intent(inout) :: w
        type(arena_t), intent(inout), target :: arena
        character(*), intent(in) :: generator, path

        w%arena => arena
        w%engine = make_native_engine(arena)
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

    !> A binding of an opaque applied function (or one of its partial
    !> derivatives, by argument index) to a kernel argument name.
    function binding(head, indices, replacement) result(b)
        character(*), intent(in) :: head, replacement
        integer, intent(in) :: indices(:)
        type(kernel_binding_t) :: b

        b%head = str(head)
        allocate (b%derivative_indices(size(indices)))
        b%derivative_indices = indices
        b%replacement = str(replacement)
    end function

    !> Emit one pure elemental leaf subroutine without a module wrapper, so the
    !> consumer can include it next to its hot caller.
    subroutine writer_emit(w, name, args, outputs, roots, bindings)
        class(leaf_writer_t), intent(inout) :: w
        character(*), intent(in) :: name, args(:), outputs(:)
        type(expr_t), intent(in) :: roots(:)
        type(kernel_binding_t), intent(in), optional :: bindings(:)
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
        if (present(bindings)) spec%bindings = bindings
        source = emit_kernel(roots, spec, ok, message=message)
        if (.not. ok) error stop message
        write (w%unit, '(a)') chars(source)
    end subroutine
end module libneo_gen_support
