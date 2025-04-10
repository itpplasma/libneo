!> This is the main module of the scheduler class
!> Used to inhert a special scheduler which is specific for a special physical problem
module scheduler_module

  use list_module
  use intlist_module
  use mpiprovider_module
  use matrix_module
  use clientStatus_module
  use workunit_module
  use wuDataRequester_module
  use initWorkunit_module

  implicit none
  
  ! Header definitions
  include "scheduler_header.f90"
  
contains
  
  ! Generic part, in most cases there is no need to adapt this part for special problems
  include "scheduler_generic.f90"

  ! Specific part, these functions have to be overwritten to adapt the library to a special problem
  include "scheduler_specific.f90"
    
end module scheduler_module
