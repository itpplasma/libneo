!> Module for class commandline_parser
module commandline_parser_module
  
  !> Class definition of commandline_parser
  type commandline_parser
   contains
     procedure :: getInt => commandline_parser_getInt
  end type commandline_parser
  
  !Singleton!
  type(commandline_parser) :: comlineParser
  
contains
  
  !> Get integer value from commandline arguments
  function commandline_parser_getInt(this, argStr, defValue) result(res)
    class(commandline_parser) :: this
    character(len=*) :: argStr        !< Key string
    character(len=512) :: str
    integer :: defValue               !< Default value, if not defined
    integer :: res
    integer :: i, j
    
    res = defValue

    ! Read command line arguments
    do i = 1, command_argument_count()
       call get_command_argument(i, str)
       j = index(str, argStr)
       if (j > 0) then
          read (str(j+len(argStr):),*) res
          exit
       end if
    end do
    
  end function commandline_parser_getInt
  
end module commandline_parser_module
