!> Module for class generic workunit
module genericWorkunit_module

  use workunit_module
  use mpiprovider_module

  !> Class generic Workunit
  type, extends(workunit) :: genericWorkunit
  contains

    procedure :: init => init_genericWorkunit       !< Constructor
    procedure :: pack => pack_genericWorkunit       !< Generic pack
    procedure :: unpack => unpack_genericWorkunit   !< Generic unpack
    procedure :: deinit => deinit_genericWorkunit   !< Destructor
    procedure :: free => free_genericworkunit       !< Generic free
    procedure :: print => print_genericWorkunit     !< Empty print

    procedure :: process => process_genericWorkunit !< Represents the job of the workunit

    procedure :: setClient => setClient_genericworkunit
    procedure :: setOldClient => setOldClient_genericworkunit
    procedure :: setDrUID => setDrUID_genericworkunit

  end type genericWorkunit

contains

  !> This subroutine has to be inherited to give the workunit the opportunity to print itself
  subroutine print_genericWorkunit(this)
    class(genericWorkunit) :: this

  end subroutine print_genericWorkunit

  !> This subroutine has to be inherted to define the job of the workunit
  subroutine process_genericWorkunit(this)
    class(genericWorkunit) :: this

  end subroutine process_genericWorkunit

  !> Constructor, call this%genericWorkunit%init() in the constructor of the inherted workunit
  subroutine init_genericWorkunit(this)
    class(genericWorkunit) :: this

    this%uid = mpro%storage%nextUID
    mpro%storage%nextUId = mpro%storage%nextUID + 1

    call this%neededWUs%init()
  end subroutine init_genericWorkunit

  !> Generic pack of internal variables of the workunit. 
  !> This has to be called in the inherited workunit before packing other variables
  subroutine pack_genericworkunit(this)
    class(genericworkunit) :: this

    call mpro%packBuffer%add_string(this%type)

    call mpro%packBuffer%add_int(this%uid)
    call mpro%packBuffer%add_int(this%client)
    call mpro%packBuffer%add_bool(this%sendBack)
    call mpro%packBuffer%add_bool(this%isProcessed)

  end subroutine pack_genericworkunit

  !> Generic unpack of internal variables, has to be called before unpacking other data
  subroutine unpack_genericworkunit(this)
    class(genericworkunit) :: this

    call mpro%packBuffer%resetPos()
    call mpro%packBuffer%get_string(this%type)
    call mpro%packBuffer%get_int(this%uid)

    call mpro%packBuffer%get_int(this%client)

    call mpro%packBuffer%get_bool(this%sendBack)
    call mpro%packBuffer%get_bool(this%isProcessed)
  end subroutine unpack_genericworkunit

  !> Frees the neededWUs-list
  subroutine free_genericworkunit(this)
    class(genericworkunit) :: this

    call this%neededWUs%free()
  end subroutine free_genericworkunit

  !> Destructor
  subroutine deinit_genericworkunit(this)
    class(genericworkunit) :: this

    call this%neededWUs%free()
  end subroutine deinit_genericworkunit

  !> Setter for datarequester-id
  subroutine setDrUID_genericworkunit(this, druid)
    class(genericworkunit) :: this
    integer :: druid

    this%druid = druid
  end subroutine setDrUID_genericworkunit

  !> Setter for oldClient, used for dataRequester
  subroutine setOldClient_genericworkunit(this, cl)
    class(genericworkunit) :: this
    integer :: cl

    this%oldClient = cl
  end subroutine setOldClient_genericworkunit

  !> Setter of client id
  subroutine setClient_genericworkunit(this, cl)
    class(genericworkunit) :: this
    integer :: cl

    this%client=  cl
  end subroutine setClient_genericworkunit

end module genericWorkunit_module
