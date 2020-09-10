!> Module for class configFile_parser
module configFile_parser_module
    implicit none

    integer, parameter :: maxCount = 100

    !> Class definition of configFileParser
    type :: configFileParser
      character(len=256) :: filename      !< Configfile name
      integer :: fileHandle = 50          !< Fortran file handle
      integer :: count                    !< Count of lines
      character(80), dimension(maxCount) :: keysAndValues   !< Storage for keys and values in the file

      contains
        procedure :: readOut => readOut_configFileParser
        procedure :: getInt  => getInt_configFileParser
        procedure, private :: searchKey => searchKey_configFileParser

    end type configFileParser

    ! Singleton !
    type(configFileParser) :: cfp

  contains

    !> Read out all keys and values from the config file
    subroutine readOut_configFileParser(this)
      class(configFileParser) :: this
      integer :: f, stat
      
      write (*,*) "CONFIG FILE PARSER DEPRECATED! STOPPING PROGRAM"
      stop
      this%fileName = "./config.txt"
      f = this%fileHandle
      open(f, file=this%filename, action='read', iostat = stat)
      !write (*,*) f, stat

      if (stat == 0) then

        this%count = 0
        do
          read(f, '(A)', iostat=stat) this%keysAndValues(this%count + 1)
          if (stat == 0) then
            this%count = this%count + 1
          else

            exit
          end if
        end do

        if (stat > 0) then
          write (*,*) "An error occured, while reading the config file", stat
          stop
        end if

        close(f)

      end if
    end subroutine readOut_configFileParser

    !> Internal function to search a key
    function searchKey_configFileParser(this, key, found) result(res)
      class(configFileParser) :: this
      character(len=*) :: key         !< Key string to search
      character(len=255) :: res
      logical :: found                !< Indicator, if found
      integer :: i

      res = ""
      found = .false.

      do i = 1, maxCount
        if (index(trim(this%keysAndValues(i)), key // "=") == 1) then
          res = this%keysAndValues(i)
          found = .true.
          exit
        end if
      end do

    end function searchKey_configFileParser

    !> Get integer value out of config file
    function getInt_configFileParser(this, key, defaultVal) result(res)
      class(configFileParser) :: this
      character(len=*) :: key           !< Key string
      character(len=255) :: value
      integer :: res
      integer :: defaultVal        !< Default value
      integer :: i
      logical :: found

      res = defaultVal

      value = this%searchKey(key ,found)

      if (found) then
        read (value(len(key)+2:),*) res
        !write (*,*) value(len(key)+2:)
        !write (*,*) res
      end if

    end function getInt_configFileParser

end module configFile_parser_module
