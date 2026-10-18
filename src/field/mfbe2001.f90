!> Reader and writer for MFBE2001 cylindrical-grid field files (HINT2 output).
!>
!> Fortran sequential unformatted, big-endian, 4-byte record markers.
!> Record 1: enf, enr, enz, ena, enp, enfd (6 float64 = nplanes, nr, nz,
!> ncomp = 3, nfp, nplanes).  Record 2: r00, z00, dr0, dz0 (box centre and
!> half widths).  Then nplanes records, each bfb(1:3, 1:nr, 1:nz) =
!> (B_phi, B_R, B_Z).  The planes cover one field period,
!> phi_k = k * 2 pi / (nfp * nplanes), k = 0..nplanes-1.
!> Nodes are R_i = r00 - dr0 + i * 2 dr0 / nr, Z_j = z00 - dz0 + j * 2 dz0 / nz
!> (i, j = 0..n-1).  B = 0 marks nodes outside the computational domain.
module neo_mfbe2001
    use neo_cyl_grid_field, only: cyl_grid_field_t
    implicit none
    private
    public :: read_mfbe2001, write_mfbe2001

    integer, parameter :: dp = kind(1.0d0)

contains

    subroutine read_mfbe2001(filename, field)
        character(*), intent(in) :: filename
        type(cyl_grid_field_t), intent(out) :: field
        real(dp) :: h(6), g(4)
        real(dp), allocatable :: r(:), z(:), b(:, :, :, :), plane(:, :, :)
        integer :: unit, ios, nf, nr, nz, nfp, i, k

        open (newunit=unit, file=filename, form='unformatted', &
              access='sequential', convert='big_endian', action='read', &
              status='old', iostat=ios)
        if (ios /= 0) error stop 'mfbe2001: cannot open file'
        read (unit, iostat=ios) h
        if (ios /= 0) error stop 'mfbe2001: cannot read header record'
        read (unit, iostat=ios) g
        if (ios /= 0) error stop 'mfbe2001: cannot read box record'
        nf = nint(h(1))
        nr = nint(h(2))
        nz = nint(h(3))
        nfp = nint(h(5))
        if (nint(h(4)) /= 3 .or. nint(h(6)) /= nf) &
            error stop 'mfbe2001: unexpected header'
        allocate (r(nr), z(nz), b(3, nr, nz, nf), plane(3, nr, nz))
        do i = 1, nr
            r(i) = g(1) - g(3) + real(i - 1, dp)*2.0_dp*g(3)/real(nr, dp)
        end do
        do i = 1, nz
            z(i) = g(2) - g(4) + real(i - 1, dp)*2.0_dp*g(4)/real(nz, dp)
        end do
        do k = 1, nf
            read (unit, iostat=ios) plane
            if (ios /= 0) error stop 'mfbe2001: truncated file'
            b(1, :, :, k) = plane(2, :, :)
            b(2, :, :, k) = plane(1, :, :)
            b(3, :, :, k) = plane(3, :, :)
        end do
        close (unit)
        call field%cyl_grid_field_init(r, z, nfp, b)
    end subroutine read_mfbe2001

    subroutine write_mfbe2001(filename, field)
        character(*), intent(in) :: filename
        type(cyl_grid_field_t), intent(in) :: field
        real(dp), allocatable :: plane(:, :, :)
        integer :: unit, k

        open (newunit=unit, file=filename, form='unformatted', &
              access='sequential', convert='big_endian', action='write', &
              status='replace')
        write (unit) real(field%nphi, dp), real(field%nr, dp), &
            real(field%nz, dp), 3.0_dp, real(field%nfp, dp), &
            real(field%nphi, dp)
        ! box centre and half width reproduce R_i = r00 - dr0 + i*2*dr0/nr
        write (unit) field%r_min + 0.5_dp*field%dr*field%nr, &
            field%z_min + 0.5_dp*field%dz*field%nz, &
            0.5_dp*field%dr*field%nr, 0.5_dp*field%dz*field%nz
        allocate (plane(3, field%nr, field%nz))
        do k = 1, field%nphi
            plane(1, :, :) = field%b(2, :, :, k)
            plane(2, :, :) = field%b(1, :, :, k)
            plane(3, :, :) = field%b(3, :, :, k)
            write (unit) plane
        end do
        close (unit)
    end subroutine write_mfbe2001

end module neo_mfbe2001
