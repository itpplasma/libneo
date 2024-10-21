module neo_mesh_2d
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    
    type :: mesh_2d_t
        real(dp), dimension(:), allocatable :: x1, x2
        real(dp), dimension(:,:), allocatable :: value
        integer :: n1, n2
        real(dp) :: dx1, dx2
        logical :: is_periodic(2)
        contains
        procedure :: mesh_init
        procedure :: mesh_deinit
    end type mesh_2d_t
    
    contains
    
    subroutine mesh_init(self, x1, x2, value, is_periodic)
    
        class(mesh_2d_t), intent(out) :: self
        real(dp), dimension(:), intent(in) :: x1, x2
        real(dp), dimension(:,:), intent(in), optional :: value
        logical, intent(in), optional :: is_periodic(2)
    
        integer :: n1, n2
    
        n1 = size(x1)
        n2 = size(x2)
        if (allocated(self%x1)) deallocate(self%x1)
        if (allocated(self%x2)) deallocate(self%x2)
        if (allocated(self%value)) deallocate(self%value)
    
        allocate(self%x1(n1), self%x2(n2))
        allocate(self%value(n1,n2))
        self%x1 = x1
        self%x2 = x2
        self%n1 = n1
        self%n2 = n2
        self%dx1 = self%x1(2) - self%x1(1)
        self%dx2 = self%x2(2) - self%x2(1)
        if (present(value)) then
            self%value = value
        else
            self%value = 0.0_dp
        end if
        if (present(is_periodic)) then
            self%is_periodic = is_periodic
        else
            self%is_periodic = .false.
        end if
    end subroutine mesh_init
    
    subroutine mesh_deinit(self)
        class(mesh_2d_t), intent(inout) :: self
    
        deallocate(self%x1, self%x2)
        deallocate(self%value)
        self%n1 = 0
        self%n2 = 0
        self%dx1 = 0.0_dp
        self%dx2 = 0.0_dp
        self%is_periodic = .false.
    end subroutine mesh_deinit
        
    end module neo_mesh_2d