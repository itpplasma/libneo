program h5merge

  !************************************
  ! HDF5
  !************************************
  USE hdf5_tools
  USE hdf5_tools_f2003

  implicit none

  integer :: funit = 10
  integer :: stat
  character(256) :: surfname, collectionname, cwd

  integer(HID_T) :: h5id_final, h5id_collection, h5id_group

  call h5_init()

  call getcwd(cwd)
  write (collectionname,*) trim(adjustl(cwd(index(cwd, '/', .true.)+1:)))
  call h5_open_rw(trim(adjustl(collectionname)) // '.h5', h5id_collection)

  open(funit, file='jobs_list.txt', status='old', action='read')
  do
     read(funit, *, iostat=stat) surfname
     if (stat /= 0) exit

     write (*,*) "Processing ", trim(surfname)
     call h5_open(trim(surfname) // '/final.h5', h5id_final)
     call h5_open_group(h5id_final, trim(surfname), h5id_group)
     !call h5_define_group(h5id_collection, trim(surfname), h5id_group)
     call h5_copy(h5id_group, '.' , h5id_collection, '/' // trim(surfname))
     call h5_close_group(h5id_group)
     call h5_close(h5id_final)
  end do

  call h5_close(h5id_collection)

  close(funit)

  call h5_deinit()

end program h5merge
