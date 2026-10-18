program test_eqfile_empty_boundaries
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_sub, only: read_eqfile2
    use input_files, only: gfile
    implicit none

    integer :: nr, nz, unit, i, j, count_boundary
    real(dp) :: psi_axis, psi_edge, b0, r0, f(4), r(4), z(4), psi(4, 4)

    do count_boundary = 0, 4, 4
        gfile = 'empty_boundary.geqdsk'
        open (newunit=unit, file=trim(gfile), status='replace')
        write (unit, '(a48,3i4)') 'No limiter or trailing empty record', 0, 4, 4
        write (unit, '(5e16.9)') 2.0_dp, 2.0_dp, 10.0_dp, 9.0_dp, 0.0_dp
        write (unit, '(5e16.9)') 10.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 2.0_dp
        write (unit, '(5e16.9)') 1.0_dp, 0.0_dp, 0.0_dp, 10.0_dp, 0.0_dp
        write (unit, '(5e16.9)') 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp
        write (unit, '(5e16.9)') (20.0_dp, i=1, 4)
        do j = 1, 3
            write (unit, '(5e16.9)') (0.0_dp, i=1, 4)
        end do
        write (unit, '(5e16.9)') ((real(10*i+j, dp), i=1, 4), j=1, 4)
        write (unit, '(5e16.9)') (2.0_dp, i=1, 4)
        write (unit, '(2i5)') count_boundary, 0
        if (count_boundary > 0) then
            write (unit, '(5e16.9)') (real(i, dp), i=1, 2*count_boundary)
        end if
        close (unit)
        nr = 4
        nz = 4
        call read_eqfile2(nr, nz, psi_axis, psi_edge, b0, r0, f, r, z, psi)
        if (abs(psi(2, 3)-23.0_dp) > 1e-12_dp) error stop 'Incorrect flux'
        if (any(abs(f-20.0_dp) > 1e-12_dp)) error stop 'Incorrect toroidal field'
        if (abs(r(1)-9.0_dp) > 1e-12_dp) error stop 'Incorrect R grid'
        open (newunit=unit, file=trim(gfile), status='old')
        close (unit, status='delete')
    end do
    print *, 'empty-boundary reader returned correct fields'
end program test_eqfile_empty_boundaries
