submodule(libneo_coordinates) libneo_coordinates_file_detection
    use netcdf, only: NF90_NOERR, NF90_NOWRITE, nf90_close, nf90_inq_dimid, &
                      nf90_inq_varid, nf90_open, nf90_strerror
    implicit none

contains

    module subroutine detect_refcoords_file_type(filename, file_type, ierr, message)
        character(len=*), intent(in) :: filename
        integer, intent(out) :: file_type
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: message

        integer :: ncid
        integer :: status
        integer :: dim_rho, dim_theta, dim_zeta
        integer :: var_x, var_y, var_z
        integer :: var_rmnc
        integer :: ierr_val
        character(len=2048) :: message_val

        ierr = 0
        message = ""
        file_type = refcoords_file_unknown

        status = nf90_open(trim(filename), NF90_NOWRITE, ncid)
        if (status /= NF90_NOERR) then
            ierr = 1
            message = "open failed: "//trim(nf90_strerror(status))
            return
        end if

        if (nf90_inq_dimid(ncid, "rho", dim_rho) == NF90_NOERR .and. &
            nf90_inq_dimid(ncid, "theta", dim_theta) == NF90_NOERR .and. &
            nf90_inq_dimid(ncid, "zeta", dim_zeta) == NF90_NOERR .and. &
            nf90_inq_varid(ncid, "x", var_x) == NF90_NOERR .and. &
            nf90_inq_varid(ncid, "y", var_y) == NF90_NOERR .and. &
            nf90_inq_varid(ncid, "z", var_z) == NF90_NOERR) then
            status = nf90_close(ncid)
            file_type = refcoords_file_chartmap

            call validate_chartmap_file(filename, ierr_val, message_val)
            if (ierr_val /= 0) then
                ierr = ierr_val
                message = trim(message_val)
            end if
            return
        end if

        if (nf90_inq_varid(ncid, "rmnc", var_rmnc) == NF90_NOERR) then
            file_type = refcoords_file_vmec_wout
            status = nf90_close(ncid)
            return
        end if

        status = nf90_close(ncid)
    end subroutine detect_refcoords_file_type

end submodule libneo_coordinates_file_detection
