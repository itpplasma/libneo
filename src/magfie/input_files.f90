module input_files
  character(len=1024) :: eqfile, cfile, gfile,pfile,convexfile,fluxdatapath
  integer :: iunit=1738
  integer :: ieqfile=1 ! equilibrium file type (0 - old, 1 - GEQDSK, 2 - WEST GEQDSK)

  data eqfile  /'d3d.equ'/
  data cfile   /'DATA/ccoil.dat'/
end module input_files
