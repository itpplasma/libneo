module packBuffer_module

  use mpi
  use mpelog_module

  implicit none

  INTEGER, PARAMETER    :: longint = 8
  integer, parameter :: maxStrLen = 1024


  !> Class for send- or receive-buffer
  type :: packBuffer
  
    !> Default send- and receive-buffer size
    integer :: buffersize

    integer, allocatable, dimension(:) :: sendBuffer, recvBuffer
    integer :: sendPos, recvPos
  contains

    !> Constructor of class buffer
    !procedure, nopass :: init => init_packBuffer

    !> Subroutines to write in the buffer
    procedure :: add_mat2dim_int => add_mat2dim_int_packBuffer
    procedure :: add_mat2dim_longint => add_mat2dim_longint_packBuffer
    procedure :: add_mat2dim_real => add_mat2dim_real_packBuffer
    procedure :: add_mat2dim_double => add_mat2dim_double_packBuffer
    procedure :: add_int   => add_int_packBuffer
    procedure :: add_bool  => add_bool_packBuffer
    procedure :: add_real  => add_real_packBuffer
    procedure :: add_double => add_double_packBuffer
    procedure :: add_string => add_string_packBuffer
    procedure :: add_array_int    => add_array_int_packBuffer
    procedure :: add_array_longint => add_array_longint_packBuffer
    procedure :: add_array_real   => add_array_real_packBuffer
    procedure :: add_array_double => add_array_double_packBuffer

    !> Generic add for pack buffer
    generic, public  :: add => add_mat2dim_longint, add_mat2dim_real, add_mat2dim_double, &
                               add_int, add_bool, add_real, add_double, add_string, &
                               add_array_int, add_array_double, add_array_real, add_array_longint, &
                               add_mat2dim_int

    !> Functions to read from the buffer
    procedure :: get_int   => get_int_packBuffer
    procedure :: get_bool  => get_bool_packBuffer
    procedure :: get_mat2dim_real => get_mat2dim_real_packBuffer
    procedure :: get_mat2dim_double => get_mat2dim_double_packBuffer
    procedure :: get_mat2dim_int => get_mat2dim_int_packBuffer
    procedure :: get_mat2dim_longint => get_mat2dim_longint_packBuffer
    procedure :: get_real  => get_real_packBuffer
    procedure :: get_double => get_double_packBuffer
    procedure :: get_string => get_string_packBuffer
    procedure :: get_array_int => get_array_int_packBuffer
    procedure :: get_array_longint => get_array_longint_packBuffer
    procedure :: get_array_real => get_array_real_packBuffer
    procedure :: get_array_double => get_array_double_packBuffer

    generic, public :: get => get_int, get_bool, get_mat2dim_real, get_mat2dim_double, get_mat2dim_longint, &
                              get_real, get_double, get_string, get_array_int, get_array_longint, &
                              get_array_real, get_array_double, get_mat2dim_int

    !> Clear the send buffer
    procedure :: clear => clear_packBuffer

    !> Reset reading position of receive buffer
    procedure :: resetPos => resetPos_packBuffer

    !> Wrapper functions for the MPI commands
    procedure :: sendTo  => sendTo_packBuffer
    procedure :: isendTo => isendTo_packBuffer
    procedure :: ssendTo => ssendTo_packBuffer
    procedure :: receiveFrom => receiveFrom_packBuffer

    procedure :: allocateBuffers => allocateBuffers_packBuffer
  end type packBuffer

  interface packBuffer
    procedure init_packBuffer
  end interface packBuffer

  contains

    !> Constructor of buffer
  function init_packBuffer() result(obj)
    type(packBuffer) :: obj

    call obj%resetPos()
  end function init_packBuffer

  subroutine allocateBuffers_packBuffer(this, buffersize)
    class(packBuffer) :: this
    integer :: buffersize

    this%buffersize = buffersize * 1024 * 1024 / 4

    allocate(this%sendBuffer(this%buffersize))
    allocate(this%recvBuffer(this%buffersize))

  end subroutine allocateBuffers_packBuffer

  !> Reset read position from receive buffer
  subroutine resetPos_packBuffer(this)
    class(packBuffer):: this

    this%sendPos = 0
    this%recvPos = 0
  end subroutine resetPos_packBuffer

  !> Clear send buffer
  subroutine clear_packBuffer(this)
    class(packBuffer) :: this

    this%sendPos = 0
  end subroutine clear_packBuffer

  !> Add array to buffer
  subroutine add_array_int_packBuffer(this, val)
    implicit none
    class(packBuffer) :: this
    integer, dimension(:), allocatable :: val
    integer :: ierr, ub, lb

    call MPI_Pack(allocated(val), 1, MPI_LOGICAL, this%sendBuffer, this%buffersize, this%sendPos, MPI_COMM_WORLD, ierr)
    if (allocated(val)) then
      ub = ubound(val, 1)
      lb = lbound(val, 1)

      call MPI_Pack(lb, 1, MPI_INTEGER, this%sendBuffer, this%buffersize, this%sendPos, MPI_COMM_WORLD, ierr)
      call MPI_Pack(ub, 1, MPI_INTEGER, this%sendBuffer, this%buffersize, this%sendPos, MPI_COMM_WORLD, ierr)
      call MPI_Pack(val, size(val), MPI_INTEGER, this%sendBuffer, this%buffersize, this%sendPos, MPI_COMM_WORLD, ierr)
    end if
  end subroutine add_array_int_packBuffer

  !> Add array to buffer
  subroutine add_array_longint_packBuffer(this, val)
    implicit none
    class(packBuffer) :: this
    integer(kind=longint), dimension(:), allocatable :: val
    integer :: ierr, ub, lb

    call MPI_Pack(allocated(val), 1, MPI_LOGICAL, this%sendBuffer, this%buffersize, this%sendPos, MPI_COMM_WORLD, ierr)
    if (allocated(val)) then
      ub = ubound(val, 1)
      lb = lbound(val, 1)

      call MPI_Pack(lb, 1, MPI_INTEGER, this%sendBuffer, this%buffersize, this%sendPos, MPI_COMM_WORLD, ierr)
      call MPI_Pack(ub, 1, MPI_INTEGER, this%sendBuffer, this%buffersize, this%sendPos, MPI_COMM_WORLD, ierr)
      call MPI_Pack(val, size(val), MPI_INTEGER8, this%sendBuffer, this%buffersize, this%sendPos, MPI_COMM_WORLD, ierr)
    end if
  end subroutine add_array_longint_packBuffer

  !> Add array to buffer
  subroutine add_array_double_packBuffer(this, val)
    implicit none
    class(packBuffer) :: this
    double precision, dimension(:), allocatable :: val
    integer :: ierr, ub, lb

    call MPI_Pack(allocated(val), 1, MPI_LOGICAL, this%sendBuffer, this%buffersize, this%sendPos, MPI_COMM_WORLD, ierr)
    if (allocated(val)) then
      ub = ubound(val, 1)
      lb = lbound(val, 1)

      call MPI_Pack(lb, 1, MPI_INTEGER, this%sendBuffer, this%buffersize, this%sendPos, MPI_COMM_WORLD, ierr)
      call MPI_Pack(ub, 1, MPI_INTEGER, this%sendBuffer, this%buffersize, this%sendPos, MPI_COMM_WORLD, ierr)
      call MPI_Pack(val, size(val), MPI_DOUBLE_PRECISION, this%sendBuffer, this%buffersize, this%sendPos, MPI_COMM_WORLD, ierr)
    end if
  end subroutine add_array_double_packBuffer

  !> Add array to buffer
  subroutine add_array_real_packBuffer(this, val)
    implicit none
    class(packBuffer) :: this
    real, dimension(:), allocatable :: val
    integer :: ierr, ub, lb

    call MPI_Pack(allocated(val), 1, MPI_LOGICAL, this%sendBuffer, this%buffersize, this%sendPos, MPI_COMM_WORLD, ierr)
    if (allocated(val)) then
      ub = ubound(val, 1)
      lb = lbound(val, 1)

      call MPI_Pack(lb, 1, MPI_INTEGER, this%sendBuffer, this%buffersize, this%sendPos, MPI_COMM_WORLD, ierr)
      call MPI_Pack(ub, 1, MPI_INTEGER, this%sendBuffer, this%buffersize, this%sendPos, MPI_COMM_WORLD, ierr)
      call MPI_Pack(val, size(val), MPI_REAL, this%sendBuffer, this%buffersize, this%sendPos, MPI_COMM_WORLD, ierr)
    end if
  end subroutine add_array_real_packBuffer


  !> Add 2dim - matrix to buffer
  subroutine add_mat2dim_real_packBuffer(this, val)
    class(packBuffer) :: this
    real, dimension(:,:), allocatable :: val
    integer :: ierr, ub1, ub2, lb1, lb2

    call MPI_Pack(allocated(val), 1, MPI_LOGICAL, this%sendBuffer, this%buffersize, this%sendPos, MPI_COMM_WORLD, ierr)
    if (allocated(val)) then
      lb1 = lbound(val, 1)
      ub1 = ubound(val, 1)

      lb2 = lbound(val, 2)
      ub2 = ubound(val, 2)

      call MPI_Pack(lb1, 1, MPI_INTEGER, this%sendBuffer, this%buffersize, this%sendPos, MPI_COMM_WORLD, ierr)
      call MPI_Pack(ub1, 1, MPI_INTEGER, this%sendBuffer, this%buffersize, this%sendPos, MPI_COMM_WORLD, ierr)
      call MPI_Pack(lb2, 1, MPI_INTEGER, this%sendBuffer, this%buffersize, this%sendPos, MPI_COMM_WORLD, ierr)
      call MPI_Pack(ub2, 1, MPI_INTEGER, this%sendBuffer, this%buffersize, this%sendPos, MPI_COMM_WORLD, ierr)

      call MPI_PACK(val, size(val), MPI_REAL, this%sendBuffer, this%buffersize, this%sendPos, MPI_COMM_WORLD, ierr)
    end if
  end subroutine add_mat2dim_real_packBuffer

  !> Add 2dim - matrix to buffer
  subroutine add_mat2dim_double_packBuffer(this, val)
    class(packBuffer) :: this
    double precision, dimension(:,:), allocatable :: val
    integer :: ierr, ub1, ub2, lb1, lb2

    call MPI_Pack(allocated(val), 1, MPI_LOGICAL, this%sendBuffer, this%buffersize, this%sendPos, MPI_COMM_WORLD, ierr)
    if (allocated(val)) then
      lb1 = lbound(val, 1)
      ub1 = ubound(val, 1)

      lb2 = lbound(val, 2)
      ub2 = ubound(val, 2)

      call MPI_Pack(lb1, 1, MPI_INTEGER, this%sendBuffer, this%buffersize, this%sendPos, MPI_COMM_WORLD, ierr)
      call MPI_Pack(ub1, 1, MPI_INTEGER, this%sendBuffer, this%buffersize, this%sendPos, MPI_COMM_WORLD, ierr)
      call MPI_Pack(lb2, 1, MPI_INTEGER, this%sendBuffer, this%buffersize, this%sendPos, MPI_COMM_WORLD, ierr)
      call MPI_Pack(ub2, 1, MPI_INTEGER, this%sendBuffer, this%buffersize, this%sendPos, MPI_COMM_WORLD, ierr)

      call MPI_PACK(val, size(val), MPI_DOUBLE_PRECISION, this%sendBuffer, this%buffersize, this%sendPos, MPI_COMM_WORLD, ierr)
    end if

  end subroutine add_mat2dim_double_packBuffer

  subroutine add_mat2dim_int_packBuffer(this, val)
    class(packBuffer) :: this
    integer, dimension(:,:), allocatable :: val
    integer :: ub1, ub2, lb1, lb2, ierr

    call MPI_Pack(allocated(val), 1, MPI_LOGICAL, this%sendBuffer, this%buffersize, this%sendPos, MPI_COMM_WORLD, ierr)
    if (allocated(val)) then
       lb1 = lbound(val, 1)
       ub1 = ubound(val, 1)

       lb2 = lbound(val, 2)
       ub2 = ubound(val, 2)

       call MPI_Pack(lb1, 1, MPI_INTEGER, this%sendBuffer, this%buffersize, this%sendPos, MPI_COMM_WORLD, ierr)
       call MPI_Pack(ub1, 1, MPI_INTEGER, this%sendBuffer, this%buffersize, this%sendPos, MPI_COMM_WORLD, ierr)
       call MPI_Pack(lb2, 1, MPI_INTEGER, this%sendBuffer, this%buffersize, this%sendPos, MPI_COMM_WORLD, ierr)
       call MPI_Pack(ub2, 1, MPI_INTEGER, this%sendBuffer, this%buffersize, this%sendPos, MPI_COMM_WORLD, ierr)

       call MPI_PACK(val, size(val), MPI_INTEGER, this%sendBuffer, this%buffersize, this%sendPos, MPI_COMM_WORLD, ierr)
    end if
  end subroutine add_mat2dim_int_packBuffer

  subroutine add_mat2dim_longint_packBuffer(this, val)
    class(packBuffer) :: this
    integer(kind=longint), dimension(:,:), allocatable :: val
    integer :: ub1, ub2, lb1, lb2, ierr

    call MPI_Pack(allocated(val), 1, MPI_LOGICAL, this%sendBuffer, this%buffersize, this%sendPos, MPI_COMM_WORLD, ierr)
    if (allocated(val)) then
     lb1 = lbound(val, 1)
      ub1 = ubound(val, 1)

      lb2 = lbound(val, 2)
      ub2 = ubound(val, 2)

      call MPI_Pack(lb1, 1, MPI_INTEGER, this%sendBuffer, this%buffersize, this%sendPos, MPI_COMM_WORLD, ierr)
      call MPI_Pack(ub1, 1, MPI_INTEGER, this%sendBuffer, this%buffersize, this%sendPos, MPI_COMM_WORLD, ierr)
      call MPI_Pack(lb2, 1, MPI_INTEGER, this%sendBuffer, this%buffersize, this%sendPos, MPI_COMM_WORLD, ierr)
      call MPI_Pack(ub2, 1, MPI_INTEGER, this%sendBuffer, this%buffersize, this%sendPos, MPI_COMM_WORLD, ierr)

      call MPI_PACK(val, size(val), MPI_INTEGER8, this%sendBuffer, this%buffersize, this%sendPos, MPI_COMM_WORLD, ierr)
    end if
  end subroutine add_mat2dim_longint_packBuffer

  !> Add integer to buffer
  subroutine add_int_packBuffer(this, val)
    class(packBuffer) :: this
    integer :: val
    integer :: ierr

    call MPI_PACK(val, 1, MPI_INTEGER, this%sendBuffer, this%buffersize, this%sendPos, MPI_COMM_WORLD, ierr)
  end subroutine add_int_packBuffer

  !> Add boolean to buffer
  subroutine add_bool_packBuffer(this, val)
    class(packBuffer) :: this
    logical :: val
    integer :: ierr

    call MPI_PACK(val, 1, MPI_LOGICAL, this%sendBuffer, this%buffersize, this%sendPos, MPI_COMM_WORLD, ierr)
  end subroutine add_bool_packBuffer

  !> Add real to buffer
  subroutine add_real_packBuffer(this, val)
    class(packBuffer) :: this
    real :: val
    integer :: ierr

    call MPI_PACK(val, 1, MPI_REAL, this%sendBuffer, this%buffersize, this%sendPos, MPI_COMM_WORLD, ierr)
  end subroutine add_real_packBuffer

  !> Add double precision to buffer
  subroutine add_double_packBuffer(this, val)
    class(packBuffer) :: this
    double precision :: val
    integer :: ierr

    call MPI_PACK(val, 1, MPI_DOUBLE_PRECISION, this%sendBuffer, this%buffersize, this%sendPos, MPI_COMM_WORLD, ierr)
  end subroutine add_double_packBuffer

  !> Add string to buffer
  subroutine add_string_packBuffer(this, val)
    class(packBuffer) :: this
    character(len=*) :: val
    integer :: ierr

    call MPI_PACK(val, maxStrLen, MPI_CHARACTER, this%sendBuffer, this%buffersize, this%sendPos, MPI_COMM_WORLD, ierr)
    !write (*,*) "Adding string with length ", len(val), val
  end subroutine add_string_packBuffer

  !> Unpack string from buffer
  subroutine get_string_packBuffer(this, val)
    class(packBuffer) :: this
    character(len=maxStrLen) :: val
    integer :: ierr

    call MPI_UNPACK(this%recvBuffer, this%buffersize, this%recvPos, val, maxStrLen, MPI_CHARACTER, MPI_COMM_WORLD, ierr)

  end subroutine get_string_packBuffer

  !> Unpack integer from buffer
  subroutine get_int_packBuffer(this, val)
    class(packBuffer) :: this
    integer :: val, ierr

    call MPI_UNPACK(this%recvBuffer, this%buffersize, this%recvPos, val, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

  end subroutine get_int_packBuffer

  !> Unpack real from buffer
  subroutine get_real_packBuffer(this, val)
    class(packBuffer) :: this
    real :: val
    integer :: ierr

    call MPI_UNPACK(this%recvBuffer, this%buffersize, this%recvPos, val, 1, MPI_REAL, MPI_COMM_WORLD, ierr)

  end subroutine get_real_packBuffer

  !> Unpack double precision from buffer
  subroutine get_double_packBuffer(this, val)
    class(packBuffer) :: this
    double precision :: val
    integer :: ierr

    call MPI_UNPACK(this%recvBuffer, this%buffersize, this%recvPos, val, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)

  end subroutine get_double_packBuffer

  !> Unpack boolean from buffer
  subroutine get_bool_packBuffer(this, val)
    class(packBuffer) :: this
    logical :: val
    integer :: ierr

    call MPI_UNPACK(this%recvBuffer, this%buffersize, this%recvPos, val, 1, MPI_LOGICAL, MPI_COMM_WORLD, ierr)

  end subroutine get_bool_packBuffer

  !> Unpack 2dim - matrix from buffer
  subroutine get_mat2dim_real_packBuffer(this, val)
    class(packBuffer) :: this
    real, dimension(:,:), allocatable :: val
    integer :: ierr, lb1, lb2, ub1, ub2
    logical :: alloc

    call MPI_UnPack(this%recvBuffer, this%buffersize, this%recvPos, alloc, 1, MPI_LOGICAL, MPI_COMM_WORLD, ierr)
    if (alloc) then
      call MPI_UnPack(this%recvBuffer, this%buffersize, this%recvPos, lb1, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
      call MPI_UnPack(this%recvBuffer, this%buffersize, this%recvPos, ub1, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
      call MPI_UnPack(this%recvBuffer, this%buffersize, this%recvPos, lb2, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
      call MPI_UnPack(this%recvBuffer, this%buffersize, this%recvPos, ub2, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

      allocate(val(lb1:ub1, lb2:ub2))
      call MPI_UNPACK(this%recvBuffer, this%buffersize, this%recvPos, val, &
                     (ub1-lb1+1)*(ub2-lb2+1), MPI_REAL, MPI_COMM_WORLD, ierr)
    end if
  end subroutine get_mat2dim_real_packBuffer

  !> Unpack 2dim - matrix from buffer
  subroutine get_mat2dim_double_packBuffer(this, val)
    class(packBuffer) :: this
    double precision, dimension(:,:), allocatable :: val
    integer :: ierr, lb1, lb2, ub1, ub2
    logical :: alloc

    call MPI_UnPack(this%recvBuffer, this%buffersize, this%recvPos, alloc, 1, MPI_LOGICAL, MPI_COMM_WORLD, ierr)
    if (alloc) then
      call MPI_UnPack(this%recvBuffer, this%buffersize, this%recvPos, lb1, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
      call MPI_UnPack(this%recvBuffer, this%buffersize, this%recvPos, ub1, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
      call MPI_UnPack(this%recvBuffer, this%buffersize, this%recvPos, lb2, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
      call MPI_UnPack(this%recvBuffer, this%buffersize, this%recvPos, ub2, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

      allocate(val(lb1:ub1, lb2:ub2))
      call MPI_UNPACK(this%recvBuffer, this%buffersize, this%recvPos, val, (ub1-lb1+1)*(ub2-lb2+1), &
                        MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
    end if

  end subroutine get_mat2dim_double_packBuffer

  subroutine get_mat2dim_int_packBuffer(this, val)
    class(packBuffer) :: this
    integer, dimension(:,:), allocatable :: val
    integer :: ierr, lb1, lb2, ub1, ub2
    logical :: alloc

    call MPI_UnPack(this%recvBuffer, this%buffersize, this%recvPos, alloc, 1, MPI_LOGICAL, MPI_COMM_WORLD, ierr)
    if (alloc) then
       call MPI_UnPack(this%recvBuffer, this%buffersize, this%recvPos, lb1, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
       call MPI_UnPack(this%recvBuffer, this%buffersize, this%recvPos, ub1, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
       call MPI_UnPack(this%recvBuffer, this%buffersize, this%recvPos, lb2, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
       call MPI_UnPack(this%recvBuffer, this%buffersize, this%recvPos, ub2, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

       allocate(val(lb1:ub1, lb2:ub2))

       call MPI_UNPACK(this%recvBuffer, this%buffersize, this%recvPos, val, (ub1-lb1+1)*(ub2-lb2+1), &
            MPI_INTEGER, MPI_COMM_WORLD, ierr)
    end if

  end subroutine get_mat2dim_int_packBuffer
  
  subroutine get_mat2dim_longint_packBuffer(this, val)
    class(packBuffer) :: this
    integer(kind=longint), dimension(:,:), allocatable :: val
    integer :: ierr, lb1, lb2, ub1, ub2
    logical :: alloc

    call MPI_UnPack(this%recvBuffer, this%buffersize, this%recvPos, alloc, 1, MPI_LOGICAL, MPI_COMM_WORLD, ierr)
    if (alloc) then
      call MPI_UnPack(this%recvBuffer, this%buffersize, this%recvPos, lb1, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
      call MPI_UnPack(this%recvBuffer, this%buffersize, this%recvPos, ub1, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
      call MPI_UnPack(this%recvBuffer, this%buffersize, this%recvPos, lb2, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
      call MPI_UnPack(this%recvBuffer, this%buffersize, this%recvPos, ub2, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

      allocate(val(lb1:ub1, lb2:ub2))

      call MPI_UNPACK(this%recvBuffer, this%buffersize, this%recvPos, val, (ub1-lb1+1)*(ub2-lb2+1), &
                      MPI_INTEGER8, MPI_COMM_WORLD, ierr)
    end if

  end subroutine get_mat2dim_longint_packBuffer

  subroutine get_array_int_packBuffer(this, val)
    class(packBuffer) :: this
    integer, dimension(:), allocatable :: val
    integer :: lb, ub, ierr
    logical :: alloc

    call MPI_UnPack(this%recvBuffer, this%buffersize, this%recvPos, alloc, 1, MPI_LOGICAL, MPI_COMM_WORLD, ierr)
    if (alloc) then
      call MPI_UnPack(this%recvBuffer, this%buffersize, this%recvPos, lb, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
      call MPI_UnPack(this%recvBuffer, this%buffersize, this%recvPos, ub, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

      allocate(val(lb:ub))
      call MPI_UNPACK(this%recvBuffer, this%buffersize, this%recvPos, val, ub-lb+1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    end if
  end subroutine get_array_int_packBuffer

  subroutine get_array_longint_packBuffer(this, val)
    class(packBuffer) :: this
    integer(kind=longint), dimension(:), allocatable :: val
    integer :: lb, ub, ierr
    logical :: alloc

    call MPI_UnPack(this%recvBuffer, this%buffersize, this%recvPos, alloc, 1, MPI_LOGICAL, MPI_COMM_WORLD, ierr)
    if (alloc) then
      call MPI_UnPack(this%recvBuffer, this%buffersize, this%recvPos, lb, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
      call MPI_UnPack(this%recvBuffer, this%buffersize, this%recvPos, ub, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

      allocate(val(lb:ub))
      call MPI_UNPACK(this%recvBuffer, this%buffersize, this%recvPos, val, ub-lb+1, MPI_INTEGER8, MPI_COMM_WORLD, ierr)
    end if

  end subroutine get_array_longint_packBuffer

  subroutine get_array_real_packBuffer(this, val)
    class(packBuffer) :: this
    real, dimension(:), allocatable :: val
    integer :: lb, ub, ierr
    logical :: alloc

    call MPI_UnPack(this%recvBuffer, this%buffersize, this%recvPos, alloc, 1, MPI_LOGICAL, MPI_COMM_WORLD, ierr)
    if (alloc) then
      call MPI_UnPack(this%recvBuffer, this%buffersize, this%recvPos, lb, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
      call MPI_UnPack(this%recvBuffer, this%buffersize, this%recvPos, ub, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

      allocate(val(lb:ub))
      call MPI_UNPACK(this%recvBuffer, this%buffersize, this%recvPos, val, ub-lb+1, MPI_REAL, MPI_COMM_WORLD, ierr)
    end if

  end subroutine get_array_real_packBuffer

  subroutine get_array_double_packBuffer(this, val)
    class(packBuffer) :: this
    double precision, dimension(:), allocatable :: val
    integer :: lb, ub, ierr
    logical :: alloc

    call MPI_UnPack(this%recvBuffer, this%buffersize, this%recvPos, alloc, 1, MPI_LOGICAL, MPI_COMM_WORLD, ierr)
    if (alloc) then
      call MPI_UnPack(this%recvBuffer, this%buffersize, this%recvPos, lb, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
      call MPI_UnPack(this%recvBuffer, this%buffersize, this%recvPos, ub, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

      allocate(val(lb:ub))
      call MPI_UNPACK(this%recvBuffer, this%buffersize, this%recvPos, val, ub-lb+1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
    end if

  end subroutine get_array_double_packBuffer

  !> Warpper for sending and logging an event
  subroutine sendTo_packBuffer(this, destRank, tag)
    class(packBuffer) :: this
    integer :: destRank, ierr, tag

    call mlog%logEvent(mpe_e_sendA)
    call MPI_Send(this%sendBuffer, this%sendPos, MPI_PACKED, destRank, tag, MPI_COMM_WORLD, ierr)
    call mlog%logEvent(mpe_e_sendB, destRank, tag)

    call this%clear()

    if (ierr /= 0) then
      write (*,*) "ERROR in MPI_Send at sendTo_packBuffer", ierr
      stop
    end if
  end subroutine sendTo_packBuffer

  !> Wrapper for SSend and logging
  subroutine ssendTo_packBuffer(this, destRank, tag)
    class(packBuffer) :: this
    integer :: destRank, ierr, tag

    !write (*,*) "Sending buffer to ", destRank

    call mlog%logEvent(mpe_e_sendA)
    call MPI_SSend(this%sendBuffer, this%sendPos, MPI_PACKED, destRank, tag, MPI_COMM_WORLD, ierr)
    call mlog%logEvent(mpe_e_sendB, destRank, tag)

    call this%clear()

    if (ierr /= 0) then
      write (*,*) "ERROR in MPI_SSend at ssendTo_packBuffer", ierr
      stop
    end if
  end subroutine ssendTo_packBuffer

  !> Wrapper for ISend and logging
  function isendTo_packBuffer(this, destRank, tag) result(res)
    class(packBuffer) :: this
    integer :: destRank, ierr, tag
    integer :: res

    !write (*,*) "Sending buffer to ", destRank
    call mlog%logEvent(mpe_e_sendA)
    call MPI_ISend(this%sendBuffer, this%sendPos, MPI_PACKED, destRank, tag, MPI_COMM_WORLD, res, ierr)
    call mlog%logEvent(mpe_e_sendB, destRank, tag)
  end function isendTo_packBuffer

  !> Wrapper for receiving and loggins
  subroutine receiveFrom_packBuffer(this, sourceRank, tag)
    class(packBuffer) :: this
    integer, intent(in) :: tag
    integer :: sourceRank, ierr, source
    integer, dimension(MPI_STATUS_SIZE) :: status

     !write (*,*) "Waiting for buffer from ", sourceRank

    call mlog%logEvent(mpe_e_recvA)
    call MPI_Recv(this%recvBuffer, this%buffersize, MPI_PACKED, sourceRank, tag, MPI_COMM_WORLD, status, ierr)
    source = status(MPI_SOURCE)

    call mlog%logEvent(mpe_e_recvB, source, tag)
  end subroutine receiveFrom_packBuffer

end module packBuffer_module
