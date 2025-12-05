program test_nctools_nc_get3
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use nctools_module, only: nc_open, nc_close, nc_get
    use netcdf
    implicit none

    integer, parameter :: nrho = 3, ntheta = 4, nzeta = 2
    real(dp) :: data(nrho, ntheta, nzeta)
    real(dp) :: readback(nrho, ntheta, nzeta)
    integer :: i, j, k
    integer :: ncid, dim_r, dim_t, dim_z, varid, ierr
    character(len=*), parameter :: fname = "test_nctools_nc_get3.nc"

    do k = 1, nzeta
        do j = 1, ntheta
            do i = 1, nrho
                data(i, j, k) = 100.0_dp*real(i, dp) + 10.0_dp*real(j, dp) &
                    + real(k, dp)
            end do
        end do
    end do

    ierr = nf90_create(fname, NF90_CLOBBER, ncid)
    if (ierr /= NF90_NOERR) then
        print *, "STOP: failed to create NetCDF file in test_nctools_nc_get3"
        error stop
    end if

    ierr = nf90_def_dim(ncid, "rho", nrho, dim_r)
    ierr = nf90_def_dim(ncid, "theta", ntheta, dim_t)
    ierr = nf90_def_dim(ncid, "zeta", nzeta, dim_z)
    ierr = nf90_def_var(ncid, "R", NF90_DOUBLE, [dim_r, dim_t, dim_z], varid)
    if (ierr /= NF90_NOERR) then
        print *, "STOP: failed to define variable in test_nctools_nc_get3"
        error stop
    end if

    ierr = nf90_enddef(ncid)
    if (ierr /= NF90_NOERR) then
        print *, "STOP: failed to enddef in test_nctools_nc_get3"
        error stop
    end if

    ierr = nf90_put_var(ncid, varid, data)
    if (ierr /= NF90_NOERR) then
        print *, "STOP: failed to write data in test_nctools_nc_get3"
        error stop
    end if

    ierr = nf90_close(ncid)
    if (ierr /= NF90_NOERR) then
        print *, "STOP: failed to close NetCDF file in test_nctools_nc_get3"
        error stop
    end if

    call nc_open(fname, ncid)
    call nc_get(ncid, "R", readback)
    call nc_close(ncid)

    if (maxval(abs(readback - data)) > 0.0_dp) then
        print *, "STOP: nc_get_double_3 readback mismatch in test_nctools_nc_get3"
        error stop
    end if

    print *, "test_nctools_nc_get3: PASS"

end program test_nctools_nc_get3

