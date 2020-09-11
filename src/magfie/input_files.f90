module input_files
  character*1024 :: eqfile, cfile, gfile,pfile,convexfile,fluxdatapath
  integer :: iunit=1738
!
  data eqfile  /'d3d.equ'/
  data cfile   /'DATA/ccoil.dat'/
!  data gfile   /'gfiles/shot115452/g115452.03525'/
!  data pfile   /'Conly/probe_g129_bfield.out'/
end module input_files
