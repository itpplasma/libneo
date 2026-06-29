module vmec_alloc_sub
  implicit none
contains
  subroutine new_allocate_vmec_stuff
!
  use new_vmec_stuff_mod
  use nctools_module
  use netcdf, only : nf90_open, NF90_NOWRITE
!
  integer :: ncid, status
  integer, dimension(2) :: lens
  integer :: len_axis
  character(len=1000) :: filename
!
  filename = netcdffile
  if (len_trim(filename) == 0) filename = 'wout.nc'
  status = nf90_open(trim(filename), NF90_NOWRITE, ncid)
  if(status /= nf90_noerr) then
    print *, "new_allocate_vmec_stuff: could not find VMEC NetCDF file"
    print *, trim(nf90_strerror(status))
    print *, trim(filename)
    error stop
  end if
!
  call nc_inq_dim(ncid, 'lmns', lens)
  call nc_inq_dim(ncid, 'raxis_cc', len_axis)
!
  call nc_get(ncid, 'nfp', nper)
!
  call nc_get(ncid, 'Rmajor_p', rmajor)
!
  nstrm=lens(1)
  nsurfm=lens(2)
  kpar=nsurfm-1
  ntor_axis=len_axis-1
  rmajor=rmajor*vmec_RZ_scale
!
  call nc_close(ncid)
!
  allocate(axm(nstrm),axn(nstrm),soa(0:kpar))
  allocate(aiota(0:kpar),sps(0:kpar),phi(0:kpar),s(0:kpar))
  allocate(rmnc(nstrm,0:kpar),zmnc(nstrm,0:kpar),almnc(nstrm,0:kpar))
  allocate(rmns(nstrm,0:kpar),zmns(nstrm,0:kpar),almns(nstrm,0:kpar))
  ! Persistent axis arrays: free a stale copy (re-entry with a different wout)
  ! before re-allocating; see new_deallocate_vmec_stuff for why they persist.
  if (allocated(raxis_cc)) deallocate(raxis_cc,zaxis_cs,raxis_cs,zaxis_cc)
  allocate(raxis_cc(0:ntor_axis),zaxis_cs(0:ntor_axis))
  allocate(raxis_cs(0:ntor_axis),zaxis_cc(0:ntor_axis))
!
  end subroutine new_allocate_vmec_stuff
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine new_deallocate_vmec_stuff
!
  use new_vmec_stuff_mod
!
  implicit none
!
  deallocate(axm,axn,soa,aiota,sps,phi,s,rmnc,zmnc,almnc,rmns,zmns,almns)
  ! raxis_cc/zaxis_cs intentionally persist: the Fourier harmonics are freed
  ! here right after the splines are built, but the chartmap export needs the
  ! exact axis later (splint_vmec_data runs off the splines, not the harmonics).
  ! They are refreshed on the next wout load (guarded re-allocation above).
!
  end subroutine new_deallocate_vmec_stuff
end  module vmec_alloc_sub
