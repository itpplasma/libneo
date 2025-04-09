!> Module for abstract class packable
module packable_module

  !> Abstract class packable
  !> If an object gets inherited from this class, it has the ability of packing and printing itself
  type, abstract :: packable
     integer :: uid = -1
   contains
     procedure(iprint), deferred :: print
     procedure(ipack),  deferred :: pack
     procedure(iunpack), deferred :: unpack
     procedure(ifree),  deferred :: free
  end type packable

  abstract interface

     subroutine ifree(this)
       import packable
       class(packable) :: this
     end subroutine ifree
     
     subroutine iprint(this)
       import packable
       class(packable) :: this
     end subroutine iprint

     subroutine ipack(this)
        import packable
        class(packable) :: this
      end subroutine ipack

      subroutine iunpack(this)
        import packable
        class(packable) :: this
      end subroutine iunpack
  end interface
    
end module packable_module
