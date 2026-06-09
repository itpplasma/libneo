program test_bessel_i
    ! Reference values from python3 mpmath (mp.dps = 40), mp.besseli(n, x).
    ! Values below 1e-280 are recorded as 0 and checked with an absolute
    ! tolerance, since they approach or pass the double-precision underflow
    ! threshold.

    use libneo_kinds, only: dp
    use neo_bessel_i, only: bessel_in, bessel_in_array
    use util_for_test, only: print_test, print_ok, print_fail

    implicit none

    call test_against_mpmath
    call test_negative_argument_parity
    call test_array_matches_elemental

contains

    subroutine test_against_mpmath
        integer, parameter :: nref = 128
        integer, dimension(nref), parameter :: n_ref = [ &
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
        0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, &
        1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, &
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, &
        5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, &
        5, 5, 5, 5, 10, 10, 10, 10, 10, 10, 10, 10, &
        10, 10, 10, 10, 10, 10, 10, 10, 50, 50, 50, 50, &
        50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, &
        100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, &
        100, 100, 100, 100, 200, 200, 200, 200, 200, 200, 200, 200, &
        200, 200, 200, 200, 200, 200, 200, 200]
        real(dp), dimension(nref), parameter :: x_ref = [ &
        1.0d-300, 1.0d-100, 1.0d-20, 1.0d-6, &
        1.0d-3, 0.1d0, 1.0d0, 1.9d0, &
        2.0d0, 2.1d0, 3.0d0, 10.0d0, &
        30.0d0, 100.0d0, 300.0d0, 700.0d0, &
        1.0d-300, 1.0d-100, 1.0d-20, 1.0d-6, &
        1.0d-3, 0.1d0, 1.0d0, 1.9d0, &
        2.0d0, 2.1d0, 3.0d0, 10.0d0, &
        30.0d0, 100.0d0, 300.0d0, 700.0d0, &
        1.0d-300, 1.0d-100, 1.0d-20, 1.0d-6, &
        1.0d-3, 0.1d0, 1.0d0, 1.9d0, &
        2.0d0, 2.1d0, 3.0d0, 10.0d0, &
        30.0d0, 100.0d0, 300.0d0, 700.0d0, &
        1.0d-300, 1.0d-100, 1.0d-20, 1.0d-6, &
        1.0d-3, 0.1d0, 1.0d0, 1.9d0, &
        2.0d0, 2.1d0, 3.0d0, 10.0d0, &
        30.0d0, 100.0d0, 300.0d0, 700.0d0, &
        1.0d-300, 1.0d-100, 1.0d-20, 1.0d-6, &
        1.0d-3, 0.1d0, 1.0d0, 1.9d0, &
        2.0d0, 2.1d0, 3.0d0, 10.0d0, &
        30.0d0, 100.0d0, 300.0d0, 700.0d0, &
        1.0d-300, 1.0d-100, 1.0d-20, 1.0d-6, &
        1.0d-3, 0.1d0, 1.0d0, 1.9d0, &
        2.0d0, 2.1d0, 3.0d0, 10.0d0, &
        30.0d0, 100.0d0, 300.0d0, 700.0d0, &
        1.0d-300, 1.0d-100, 1.0d-20, 1.0d-6, &
        1.0d-3, 0.1d0, 1.0d0, 1.9d0, &
        2.0d0, 2.1d0, 3.0d0, 10.0d0, &
        30.0d0, 100.0d0, 300.0d0, 700.0d0, &
        1.0d-300, 1.0d-100, 1.0d-20, 1.0d-6, &
        1.0d-3, 0.1d0, 1.0d0, 1.9d0, &
        2.0d0, 2.1d0, 3.0d0, 10.0d0, &
        30.0d0, 100.0d0, 300.0d0, 700.0d0]
        real(dp), dimension(nref), parameter :: i_ref = [ &
        1.00000000000000000d0, 1.00000000000000000d0, &
        1.00000000000000000d0, 1.00000000000025000d0, &
        1.00000025000001563d0, 1.00250156293409560d0, &
        1.26606587775200834d0, 2.12774019405388786d0, &
        2.27958530233606727d0, 2.44628312943618229d0, &
        4.88079258586502409d0, 2815.71662846625447d0, &
        781672297823.977490d0, 1.07375170713107382d+42, &
        4.47584736793505212d+128, 1.52959334767187374d+302, &
        0.0d0, 5.00000000000000000d-101, &
        5.00000000000000000d-21, 5.00000000000062500d-7, &
        0.000500000062500002604d0, 0.0500625260470926921d0, &
        0.565159103992485027d0, 1.44824437305488895d0, &
        1.59063685463732906d0, 1.74549980883610616d0, &
        3.95337021740260940d0, 2670.98830370125465d0, &
        768532038938.956999d0, 1.06836939033816248d+42, &
        4.46838138503695441d+128, 1.52850039023390069d+302, &
        0.0d0, 1.25000000000000000d-201, &
        1.25000000000000000d-41, 1.25000000000010417d-13, &
        1.25000010416666992d-7, 0.00125104199224175912d0, &
        0.135747669767038281d0, 0.603272432943478432d0, &
        0.688948447698738204d0, 0.783902359116081187d0, &
        2.24521244092995115d0, 2281.51896772600354d0, &
        730436828561.380356d0, 1.05238431932431057d+42, &
        4.44605815870147242d+128, 1.52522620369977688d+302, &
        0.0d0, 0.0d0, &
        2.60416666666666667d-104, 2.60416666666677517d-34, &
        2.60416677517361305d-19, 2.60525192989369689d-9, &
        0.000271463155956971875d0, 0.00748302334597045769d0, &
        0.00982567932313170232d0, 0.0127511786632194772d0, &
        0.0912064776615133485d0, 777.188286403259960d0, &
        512151465476.934970d0, 9.47009387303558125d+41, &
        4.29289057901400890d+128, 1.50250253219286880d+302, &
        0.0d0, 0.0d0, &
        2.69114445546737213d-210, 2.69114445546743330d-70, &
        2.69114451662974676d-40, 2.69175614292214153d-20, &
        2.75294803983687363d-10, 1.79054034918803492d-7, &
        3.01696387935068437d-7, 4.95995999130852991d-7, &
        0.0000194643934706129687d0, 21.8917061637233705d0, &
        145831809975.967124d0, 6.49897552472014780d+41, &
        3.78772592586686868d+128, 1.42407639992908212d+302, &
        0.0d0, 0.0d0, &
        0.0d0, 0.0d0, &
        2.92028573123431308d-230, 2.92042887162101046d-130, &
        2.93463530851183814d-80, 2.57507210369426141d-66, &
        3.35304282986064164d-65, 3.85280079999064773d-64, &
        2.19098806657918388d-56, 4.75689456072683991d-30, &
        0.000145901069164689465d0, 4.82195808559408067d+36, &
        6.95799128578292545d+126, 2.56345654772554551d+301, &
        0.0d0, 0.0d0, &
        0.0d0, 0.0d0, &
        0.0d0, 0.0d0, &
        8.47367400813807887d-189, 6.40084638918546085d-161, &
        1.08217147454986053d-158, 1.42451389046070307d-156, &
        4.45447034176550278d-141, 1.08234420174920164d-88, &
        3.94764200533342795d-40, 4.64153494161619911d+21, &
        2.92447368138262196d+121, 1.21764136670788957d+299, &
        0.0d0, 0.0d0, &
        0.0d0, 0.0d0, &
        0.0d0, 0.0d0, &
        0.0d0, 0.0d0, &
        0.0d0, 0.0d0, &
        0.0d0, 8.93536581738094807d-236, &
        6.39994862872884220d-140, 1.39540586010826159d-30, &
        4.07553713409152924d+100, 7.08020442754118915d+289]
        real(dp), parameter :: rel_tol = 1.0d-13
        real(dp), parameter :: abs_tol = 1.0d-280

        real(dp) :: val, err
        integer :: i

        call print_test("test_against_mpmath")

        do i = 1, nref
            val = bessel_in(n_ref(i), x_ref(i))
            if (i_ref(i) == 0.0d0) then
                err = abs(val)
                if (err > abs_tol) then
                    call print_fail
                    print *, "n = ", n_ref(i), ", x = ", x_ref(i)
                    print *, "got ", val, ", expected |I_n| <= ", abs_tol
                    stop 1
                end if
            else
                err = abs(val - i_ref(i))/abs(i_ref(i))
                if (err > rel_tol) then
                    call print_fail
                    print *, "n = ", n_ref(i), ", x = ", x_ref(i)
                    print *, "got ", val, ", expected ", i_ref(i)
                    print *, "rel err ", err
                    stop 1
                end if
            end if
        end do

        call print_ok
    end subroutine test_against_mpmath


    subroutine test_negative_argument_parity
        integer, parameter :: n_par(6) = [0, 1, 2, 5, 51, 200]
        real(dp), parameter :: x_par(5) = &
            [1.0d-6, 0.5d0, 2.0d0, 35.0d0, 650.0d0]

        real(dp) :: vpos, vneg, expected
        integer :: i, j

        call print_test("test_negative_argument_parity")

        do i = 1, size(n_par)
            do j = 1, size(x_par)
                vpos = bessel_in(n_par(i), x_par(j))
                vneg = bessel_in(n_par(i), -x_par(j))
                expected = vpos
                if (mod(n_par(i), 2) == 1) expected = -vpos
                if (vneg /= expected) then
                    call print_fail
                    print *, "n = ", n_par(i), ", x = ", -x_par(j)
                    print *, "got ", vneg, ", expected ", expected
                    stop 1
                end if
            end do
        end do

        call print_ok
    end subroutine test_negative_argument_parity


    subroutine test_array_matches_elemental
        integer, parameter :: nmax = 200
        real(dp), parameter :: x_arr(7) = &
            [1.0d-300, 1.0d-8, 0.3d0, 2.0d0, 7.5d0, 120.0d0, 700.0d0]
        real(dp), parameter :: rel_tol = 1.0d-13
        real(dp), parameter :: abs_tol = 1.0d-280

        real(dp) :: values(0:nmax), values_neg(0:nmax), scalar, err, expected
        integer :: i, n

        call print_test("test_array_matches_elemental")

        do i = 1, size(x_arr)
            call bessel_in_array(nmax, x_arr(i), values)
            do n = 0, nmax
                scalar = bessel_in(n, x_arr(i))
                if (abs(scalar) < abs_tol) then
                    err = abs(values(n))
                    if (err > abs_tol) then
                        call print_fail
                        print *, "n = ", n, ", x = ", x_arr(i)
                        print *, "array ", values(n), ", scalar ", scalar
                        stop 1
                    end if
                else
                    err = abs(values(n) - scalar)/abs(scalar)
                    if (err > rel_tol) then
                        call print_fail
                        print *, "n = ", n, ", x = ", x_arr(i)
                        print *, "array ", values(n), ", scalar ", scalar
                        stop 1
                    end if
                end if
            end do
            call bessel_in_array(nmax, -x_arr(i), values_neg)
            do n = 0, nmax
                expected = values(n)
                if (mod(n, 2) == 1) expected = -values(n)
                if (values_neg(n) /= expected) then
                    call print_fail
                    print *, "n = ", n, ", x = ", -x_arr(i)
                    print *, "got ", values_neg(n), ", expected ", expected
                    stop 1
                end if
            end do
        end do

        call print_ok
    end subroutine test_array_matches_elemental

end program test_bessel_i
