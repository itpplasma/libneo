!> Consistency check of a G-EQDSK file (module efit_sfl): Grad-Shafranov
!> residual, contour q against qpsi, force balance with the file's p', column
!> consistency. Exits with an error unless the file passes or --allow-poor is
!> given.
!>
!> Usage: efit_sfl_check.x <gfile> [edge_fraction] [--allow-poor] [--gs-map <file>]
!>   edge_fraction: psi_edge = simag + edge_fraction (sibry - simag), default 1
program efit_sfl_check
    use, intrinsic :: iso_fortran_env, only: dp => real64, output_unit
    use efit_sfl, only: efit_sfl_t, efit_sfl_config_t, efit_sfl_init, &
                        efit_sfl_diag_write, efit_sfl_gs_map
    use geqdsk_tools, only: geqdsk_t, geqdsk_read
    implicit none
    type(efit_sfl_t) :: sfl
    type(efit_sfl_config_t) :: cfg
    type(geqdsk_t) :: g
    character(len=1024) :: gfile, arg, mapfile
    real(dp) :: frac
    real(dp), allocatable :: res(:, :), mask(:, :)
    integer :: i, j, k, u

    if (command_argument_count() < 1) error stop &
        'usage: efit_sfl_check.x <gfile> [edge_fraction] [--allow-poor] [--gs-map <file>]'
    call get_command_argument(1, gfile)
    frac = 1.0_dp
    mapfile = ''
    k = 2
    do while (k <= command_argument_count())
        call get_command_argument(k, arg)
        select case (trim(arg))
        case ('--allow-poor')
            cfg%allow_poor_equilibrium = .true.
            cfg%override_reason = 'efit_sfl_check --allow-poor'
        case ('--gs-map')
            k = k + 1
            call get_command_argument(k, mapfile)
        case default
            read (arg, *) frac
        end select
        k = k + 1
    end do
    call geqdsk_read(g, trim(gfile))
    cfg%psi_edge_rel = frac*(g%sibry - g%simag)
    call efit_sfl_init(sfl, trim(gfile), cfg)
    call efit_sfl_diag_write(sfl, output_unit)
    if (len_trim(mapfile) > 0) then
        allocate (res(sfl%nw, sfl%nh), mask(sfl%nw, sfl%nh))
        call efit_sfl_gs_map(sfl, res, mask)
        open (newunit=u, file=trim(mapfile))
        write (u, '(a)') '# R [cm], Z [cm], (Delta* psi - gsrhs)/max|gsrhs|, inside psi_edge'
        do j = 1, sfl%nh
            do i = 1, sfl%nw
                write (u, '(4es24.15)') sfl%g%R_eqd(i), sfl%g%Z_eqd(j), res(i, j), mask(i, j)
            end do
        end do
        close (u)
    end if
end program efit_sfl_check
