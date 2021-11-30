module derived_scheduler_module
  use scheduler_module

  type, extends(scheduler) :: derived_scheduler
  contains
    procedure :: initMaster => derived_initMaster
  end type derived_scheduler
contains

  subroutine derived_initMaster(this)
    class(derived_scheduler) :: this

    write(*,*) 'Derived initMaster'
  end subroutine derived_initMaster
end module derived_scheduler_module
