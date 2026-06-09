program test_fft

    use libneo_kinds, only: dp
    use neo_fft, only: fft_r2c, fft_c2c, neo_fft_plan_t, neo_fft_plan_init
    use util_for_test, only: print_test, print_ok, print_fail

    implicit none

    call test_r2c_against_mpmath
    call test_c2c_against_mpmath
    call test_analytic_tones
    call test_amplitude_scaling
    call test_roundtrip_many_sizes

contains

    subroutine test_r2c_against_mpmath
        ! reference values from python3 + mpmath at 35 significant digits;
        ! inputs span amplitudes from 1e-9 to 1e8 within each vector
        real(dp), parameter :: x8(8) = [ &
            -4564015.8409601d0, -4.030419154507829d-5, 7.960529253397963d0, &
            -95159908.30484308d0, -0.00014838754364424612d0, 0.07145171036307936d0, &
            -1.9086859784222377d-5, -1.3183261641195765d-6]
        real(dp), parameter :: cre8(5) = [ &
            -99723916.11403131d0, 62724200.56807947d0, -4564023.801618654d0, &
            -71852232.2497029d0, 90595900.35283467d0]
        real(dp), parameter :: cim8(5) = [ &
            0.0d0, 67288208.54944782d0, -95159908.37625581d0, &
            67288224.4705445d0, 0.0d0]
        real(dp), parameter :: x12(12) = [ &
            91.52963744207916d0, 3.950867134363458d-7, 9.963046597803662d-8, &
            -2040453.3052673268d0, -1.7637980643966134d-6, -3638965.5588348038d0, &
            -172.34192174461606d0, 13.908242138261407d0, 0.006748904067921288d0, &
            9.069128318060123d0, -459142.6239035053d0, -404277.2338712193d0]
        real(dp), parameter :: cre12(7) = [ &
            -6542896.550043066d0, 2572002.7740866398d0, 248320.290201094d0, &
            459406.5022097326d0, 210660.7023865663d0, -3030617.6616188125d0, &
            5624449.691161931d0]
        real(dp), parameter :: cim12(7) = [ &
            0.0d0, 3260184.3205837924d0, -3899192.199150127d0, &
            1194239.8588096828d0, -3103933.834936095d0, 4055442.661412825d0, &
            0.0d0]
        real(dp), parameter :: x15(15) = [ &
            -3.0058221081511236d-6, -71052296.53207941d0, -89.33379086098736d0, &
            -4.785723585531545d-9, 9.484987951114614d-7, 8.661381147893285d-6, &
            -3097642.819267985d0, 0.25868066016844216d0, 789.8914226062925d0, &
            9.381252774159145d-7, 772.7841559329969d0, -0.006344484720926d0, &
            0.05997818911664177d0, -1.8802624519120826d-6, 2.7718203618329995d0]
        real(dp), parameter :: cre15(8) = [ &
            -74148462.92541933d0, -62404673.56877238d0, -48500143.93908889d0, &
            -22913384.773525316d0, 9933262.477598941d0, 32427767.323771864d0, &
            59989548.19663896d0, 68541855.74606395d0]
        real(dp), parameter :: cim15(8) = [ &
            0.0d0, 30721222.330977846d0, 49855213.6352524d0, &
            70521302.25349279d0, 68842382.41400659d0, 61533033.41495762d0, &
            43583406.55226456d0, 11827988.558741482d0]
        real(dp), parameter :: x16(16) = [ &
            3853.477197902846d0, 0.08785613192365327d0, -16.724118630446206d0, &
            -1.1592801773053862d0, -61.3467070243916d0, -50475.44484800381d0, &
            -8.968262786046008d-8, -5642700.771041143d0, -110.73068335149827d0, &
            9.136599358899822d-5, 0.8740378338453207d0, -5.913152736742393d-5, &
            0.09241814583043126d0, 4885094.054939416d0, -0.0009886480947732792d0, &
            7.725945490940598d-6]
        real(dp), parameter :: cre16(9) = [ &
            -804417.5911776769d0, 7105887.828188653d0, -7404778.700282219d0, &
            -2396525.786343685d0, 3697.343295207165d0, 2404479.091055771d0, &
            7412386.701889079d0, -7097984.301375722d0, 811748.8734899537d0]
        real(dp), parameter :: cim16(9) = [ &
            0.0d0, 6719314.660611088d0, -571383.7690774967d0, &
            3324365.552855433d0, -10477320.628411636d0, 3324463.5449522617d0, &
            -571415.4672616144d0, 6719166.896207237d0, 0.0d0]
        real(dp), parameter :: x20(20) = [ &
            0.8551622286887963d0, -0.0024494505503844088d0, -0.7493100495585514d0, &
            53113153.8644942d0, 8419646.225758852d0, -2983528.356436429d0, &
            7.042484703649934d-8, 4.296948500785236d-8, 2.119330998308535d-8, &
            -22.4020409335159d0, -0.013187600932107224d0, -64.39580417380786d0, &
            0.13115020991195703d0, -0.06050940917854732d0, 4.389774846104582d-6, &
            -0.8674396374949596d0, 3.261373721982823d-5, -0.03082413214222337d0, &
            0.010919589423488052d0, 0.7235719980158049d0]
        real(dp), parameter :: cre20(11) = [ &
            58549184.9330924d0, 33821025.7324934d0, -20241043.783378534d0, &
            -57325195.41592207d0, -43351184.7467449d0, 8419647.963677527d0, &
            48554814.95402114d0, 43701923.923454285d0, 6617771.325445109d0, &
            -28617397.861953992d0, -41709892.01203175d0]
        real(dp), parameter :: cim20(11) = [ &
            0.0d0, -47993488.32248719d0, -55462529.22823678d0, &
            -14447483.34375305d0, 39226728.93341624d0, 56096640.177082784d0, &
            23211609.37812711d0, -24345372.80125689d0, -45564642.662808195d0, &
            -31978370.554621574d0, 0.0d0]
        real(dp), parameter :: x97(97) = [ &
            -3931.4829316771884d0, -2131.7127291325955d0, -24532432.004680827d0, &
            22893634.712887857d0, 0.05318019004740961d0, -74889.87614878193d0, &
            -9.148238996187328d-7, -744.725021020733d0, 4040650.348700097d0, &
            -8.809149429009629d-8, -0.0003013066931598385d0, 61.080122247470236d0, &
            0.0008849186268537579d0, -0.09781498234310734d0, 56636.12847367405d0, &
            5.193556600430064d-7, 0.009461850469449447d0, -547532.5793472368d0, &
            -9.758793824580728d-9, 2.454918373061821d-10, 86.8462280994697d0, &
            35793.86106156077d0, 8.419079992842904d-7, 2937.358242872574d0, &
            5.41521726846085d0, -9.023209443722928d-7, 9.514237901779786d0, &
            -9.96098913378007d0, 7968.903718418059d0, -7.742336848426756d-5, &
            -7737337.2196830185d0, 67989.76949841318d0, -0.0021022145090277577d0, &
            -8.383934967921744d-9, 7.772368244239614d0, -7.36614642965201d-7, &
            7.715657315574553d-8, 7.076593497013135d-9, 4.4523488168850543d-7, &
            90490.28529890148d0, -3.8560199265855984d-6, 0.007445937798764248d0, &
            -788571.1717613526d0, 3363268.833607922d0, -21309589.143514995d0, &
            -7.609232262253318d-6, -86742397.1351574d0, 0.03821928451257386d0, &
            170.40727420187142d0, -46.350968913531275d0, 9.036258289909134d-7, &
            0.002484791424357662d0, 3.2089682931642915d-5, -3.0529496232671916d-9, &
            8.29037580634564d-9, 2.0424026507178983d0, -10332.35094158309d0, &
            5.897176239357518d-6, 3200560.6062649596d0, 0.24642618331892052d0, &
            -7.665275206756347d-7, -595153.5557148145d0, 38215.0496966609d0, &
            4449.441643311203d0, -0.07709162963612182d0, -0.0006515327538795375d0, &
            3.051398508935097d-9, -0.0004905384179246828d0, -9.80395016845242d0, &
            263.95791428247327d0, -8.38660546691349d-5, 53782300.63732556d0, &
            -0.0006372880521317115d0, 771786.954177811d0, -1.8091098416774542d-7, &
            -3.864875229707103d-6, 9.743902396681327d-6, -8.037380900820834d0, &
            75446.23475115185d0, 2.6924656462992514d-8, -3.457265346927587d0, &
            -646.7712185048023d0, -0.09674802018797672d0, 528152.1498488649d0, &
            3.3543137128694744d0, -50.83578012283889d0, 75.62503537103902d0, &
            1.6960154265545468d-8, -1.0525941616172441d-5, 23595312.563585956d0, &
            3.178394885436202d-8, 63314642.14900522d0, 6724443.1256441d0, &
            -756094.1825909519d0, 5.461575921435988d-7, 8402748.859533532d0, &
            9902.246408824818d0]
        real(dp), parameter :: cre97(49) = [ &
            47906103.958839804d0, 192794682.0671226d0, -74890405.66828686d0, &
            133201497.19961056d0, -39353542.04677024d0, -3459982.798757083d0, &
            -172764841.7553503d0, -25997458.11952775d0, -60272801.710882306d0, &
            -132065886.84077552d0, -57500995.81704479d0, -8027793.857559445d0, &
            69511643.07924539d0, -94069670.27369232d0, 111106712.33655274d0, &
            32902775.839597907d0, 93419918.08633047d0, -97673914.50489694d0, &
            123196729.81642407d0, -34303404.63707333d0, 29645707.438694987d0, &
            -112094648.68859617d0, 79833051.11849883d0, -53235254.60607848d0, &
            -19479577.86954854d0, -69269534.51776142d0, 61847799.09743962d0, &
            -20599284.53730681d0, -24419191.938595608d0, 8059304.087909336d0, &
            94952000.46337813d0, 47449748.8570593d0, 4454022.997223225d0, &
            140540906.80396217d0, 93738617.13003819d0, 89791478.36796531d0, &
            -43237702.49790945d0, 162268557.76139206d0, -70918169.15828794d0, &
            -7804452.377255638d0, -169534447.50654024d0, 66287179.8299248d0, &
            -198215527.5176798d0, -42020302.32644955d0, -121868249.11187822d0, &
            85170229.22190244d0, -81377433.0244201d0, 26750878.73915307d0, &
            63387304.46789339d0]
        real(dp), parameter :: cim97(49) = [ &
            0.0d0, 118373286.87937608d0, 8309088.771951072d0, &
            85239306.07177824d0, 45591064.15929764d0, 201309807.5807198d0, &
            -74977175.0966739d0, 86709896.03877687d0, -56502821.06616922d0, &
            82553140.60488227d0, -187037726.9478641d0, 14479775.611728273d0, &
            -86416306.95640709d0, 37886489.96023292d0, -116588463.6801031d0, &
            51567703.43013574d0, 60982474.80500791d0, 57861127.58451161d0, &
            24424105.197268635d0, 92124923.89288306d0, 149922114.34563372d0, &
            11839426.28602952d0, 59728716.81213793d0, 63782464.40458324d0, &
            106009167.35327199d0, -95389661.1225349d0, 44695266.583213404d0, &
            -28788490.596207306d0, 16606559.536430156d0, -176587142.7159333d0, &
            30869707.64733487d0, -86588492.49371903d0, -16792931.865199655d0, &
            -148613702.64922053d0, 74492760.9716158d0, -42013168.097151615d0, &
            25690979.045400202d0, -58660146.17282627d0, 91421417.69619673d0, &
            -6775993.556862573d0, -62200194.85925117d0, -41528108.453317516d0, &
            -52476476.728243805d0, -39515056.73033772d0, -235305639.59201506d0, &
            1526142.4464737908d0, -129964856.78024466d0, 13673624.895472456d0, &
            -175491763.5501008d0]

        call check_r2c("fft_r2c vs mpmath, n=8 (radix 4,2)", x8, cre8, cim8)
        call check_r2c("fft_r2c vs mpmath, n=12 (radix 2,3)", x12, cre12, cim12)
        call check_r2c("fft_r2c vs mpmath, n=15 (odd, radix 3,5)", x15, cre15, cim15)
        call check_r2c("fft_r2c vs mpmath, n=16 (radix 4,4)", x16, cre16, cim16)
        call check_r2c("fft_r2c vs mpmath, n=20 (radix 2,5)", x20, cre20, cim20)
        call check_r2c("fft_r2c vs mpmath, n=97 (prime, Bluestein)", x97, cre97, cim97)
    end subroutine test_r2c_against_mpmath

    subroutine check_r2c(name, x, cre, cim)
        character(*), intent(in) :: name
        real(dp), intent(in) :: x(:), cre(:), cim(:)
        ! tolerance: abs error per bin below 1e-14 * n * max|x|
        real(dp), parameter :: tol_fac = 1.0d-14
        complex(dp) :: c(size(cre)), cp(size(cre))
        type(neo_fft_plan_t) :: plan
        real(dp) :: tol
        integer :: n, k

        call print_test(name)
        n = size(x)
        tol = tol_fac*real(n, dp)*maxval(abs(x))
        call fft_r2c(x, c)
        call neo_fft_plan_init(plan, n)
        call fft_r2c(x, cp, plan)
        do k = 1, size(cre)
            if (abs(c(k) - cmplx(cre(k), cim(k), dp)) > tol) then
                call print_fail
                print *, "bin ", k - 1, ": got ", c(k), ", expected ", cre(k), cim(k)
                stop "fft_r2c deviates from mpmath reference"
            end if
            if (abs(cp(k) - c(k)) > 0.0d0) then
                call print_fail
                stop "planned and plan-free fft_r2c disagree"
            end if
        end do
        call print_ok
    end subroutine check_r2c

    subroutine test_c2c_against_mpmath
        ! reference values from python3 + mpmath at 35 significant digits
        real(dp), parameter :: zr10(10) = [ &
            -1.386542330195478d0, -2.554644485808514d0, -2.6758514055699862d0, &
            0.021448515339408836d0, 0.9641758091670756d0, 1.7051897448131283d0, &
            0.1170279835367376d0, 0.8145135816762661d0, -0.696313514804852d0, &
            -0.8407642047664727d0]
        real(dp), parameter :: zi10(10) = [ &
            -0.5983412370626064d0, 0.26987936667813806d0, -2.970046697163583d0, &
            0.5775171982904301d0, 0.6319626979006783d0, 1.3674675000886989d0, &
            1.2140790557626282d0, 0.7888711624067373d0, 1.3044452384117156d0, &
            -1.867411499239117d0]
        real(dp), parameter :: zre10(10) = [ &
            -4.531760306612687d0, -11.365954747495707d0, 1.8534073409711793d0, &
            5.811836444870357d0, 6.8698500373388445d0, -2.8232466091203197d0, &
            -4.055552573615914d0, -2.419709382644631d0, 1.4572925750068288d0, &
            -4.661586080652732d0]
        real(dp), parameter :: zim10(10) = [ &
            0.71842278607372d0, -2.5423002490488105d0, 4.221062138188764d0, &
            0.7458618267406261d0, 0.8525476647790614d0, -1.5542246703760545d0, &
            0.09868398332123995d0, 2.3566289848271302d0, -2.0450852572323237d0, &
            -8.835009577899418d0]
        ! tolerance: abs error per bin below 1e-14 * n * max|z|
        real(dp), parameter :: tol = 1.0d-14*10*3.0d0
        complex(dp) :: z(10)
        integer :: k

        call print_test("fft_c2c forward/backward vs mpmath, n=10")
        z = cmplx(zr10, zi10, dp)
        call fft_c2c(z, -1)
        do k = 1, 10
            if (abs(z(k) - cmplx(zre10(k), zim10(k), dp)) > tol) then
                call print_fail
                print *, "bin ", k - 1, ": got ", z(k), ", expected ", &
                    zre10(k), zim10(k)
                stop "fft_c2c deviates from mpmath reference"
            end if
        end do
        call fft_c2c(z, 1)
        z = z/10.0d0
        do k = 1, 10
            if (abs(z(k) - cmplx(zr10(k), zi10(k), dp)) > tol) then
                call print_fail
                stop "fft_c2c backward does not invert forward"
            end if
        end do
        call print_ok
    end subroutine test_c2c_against_mpmath

    subroutine test_analytic_tones
        call check_tones(360, 3, 170)
        call check_tones(1009, 5, 400)
        call check_tones(4096, 1, 2000)
    end subroutine test_analytic_tones

    subroutine check_tones(n, f1, f2)
        integer, intent(in) :: n, f1, f2
        ! x_j = a0 + a1 cos(2 pi f1 j/n + phi1) + a2 sin(2 pi f2 j/n) has the
        ! exact spectrum c(1) = n a0, c(f1+1) = n a1/2 exp(i phi1),
        ! c(f2+1) = -i n a2/2, zero elsewhere
        real(dp), parameter :: a0 = 0.5d0, a1 = 2.0d0, phi1 = 0.7d0, a2 = -1.5d0
        real(dp), parameter :: twopi = 8.0d0*atan(1.0d0)
        ! tolerance: abs error per bin below 1e-13 * n * max|x|
        real(dp), parameter :: tol_fac = 1.0d-13
        character(64) :: name
        real(dp) :: x(n), tol
        complex(dp) :: c(n/2 + 1), expected(n/2 + 1)
        type(neo_fft_plan_t) :: plan
        integer :: j, k

        write (name, '("fft_r2c analytic tones, n=", i0)') n
        call print_test(trim(name))
        ! reduce f*(j-1) mod n so the phase stays in [0, 2 pi); otherwise the
        ! rounding of large arguments perturbs the tone coherently by ~n*f*eps
        do j = 1, n
            x(j) = a0 + a1*cos(twopi*modulo(f1*(j - 1), n)/real(n, dp) + phi1) &
                   + a2*sin(twopi*modulo(f2*(j - 1), n)/real(n, dp))
        end do
        expected = (0.0d0, 0.0d0)
        expected(1) = n*a0
        expected(f1 + 1) = 0.5d0*n*a1*cmplx(cos(phi1), sin(phi1), dp)
        expected(f2 + 1) = cmplx(0.0d0, -0.5d0*n*a2, dp)
        tol = tol_fac*real(n, dp)*maxval(abs(x))
        call neo_fft_plan_init(plan, n)
        call fft_r2c(x, c, plan)
        do k = 1, n/2 + 1
            if (abs(c(k) - expected(k)) > tol) then
                call print_fail
                print *, "bin ", k - 1, ": got ", c(k), ", expected ", expected(k)
                stop "fft_r2c deviates from analytic tone spectrum"
            end if
        end do
        call print_ok
    end subroutine check_tones

    subroutine test_amplitude_scaling
        call check_scaling(256)
        call check_scaling(97)
    end subroutine test_amplitude_scaling

    subroutine check_scaling(n)
        integer, intent(in) :: n
        ! the DFT is linear and scaling by a power of two is exact in binary
        ! floating point, so fft(alpha x) = alpha fft(x) bitwise for alpha = 2^k;
        ! this covers amplitudes from 1e-120 to 1e120
        real(dp), parameter :: alpha_up = 2.0d0**400, alpha_down = 2.0d0**(-400)
        character(64) :: name
        real(dp) :: x(n)
        complex(dp) :: c(n/2 + 1), c_up(n/2 + 1), c_down(n/2 + 1)
        type(neo_fft_plan_t) :: plan
        integer :: k

        write (name, '("fft_r2c power-of-two amplitude scaling, n=", i0)') n
        call print_test(trim(name))
        call fill_random(x)
        call neo_fft_plan_init(plan, n)
        call fft_r2c(x, c, plan)
        call fft_r2c(alpha_up*x, c_up, plan)
        call fft_r2c(alpha_down*x, c_down, plan)
        do k = 1, n/2 + 1
            if (abs(c_up(k) - alpha_up*c(k)) > 0.0d0 &
                .or. abs(c_down(k) - alpha_down*c(k)) > 0.0d0) then
                call print_fail
                print *, "bin ", k - 1
                stop "fft_r2c amplitude scaling is not exact"
            end if
        end do
        call print_ok
    end subroutine check_scaling

    subroutine test_roundtrip_many_sizes
        integer, parameter :: sizes(37) = [ &
            1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 16, 20, 24, 30, 36, 49, &
            60, 64, 97, 100, 121, 128, 180, 210, 256, 360, 480, 512, 720, &
            960, 1009, 1024, 2048, 4093, 4096]
        integer :: i

        call print_test("fft_r2c/fft_c2c roundtrip over size battery")
        do i = 1, size(sizes)
            call check_roundtrip(sizes(i))
        end do
        call print_ok
    end subroutine test_roundtrip_many_sizes

    subroutine check_roundtrip(n)
        integer, intent(in) :: n
        ! rebuild the full spectrum from fft_r2c output via Hermitian symmetry,
        ! invert with the unnormalized backward fft_c2c divided by n;
        ! tolerance: abs error per sample below 1e-14 * n * max|x|
        real(dp), parameter :: tol_fac = 1.0d-14
        real(dp) :: x(n), tol
        complex(dp) :: c(n/2 + 1), z(n)
        integer :: k

        call fill_random(x)
        call fft_r2c(x, c)
        do k = 0, n/2
            z(k + 1) = c(k + 1)
        end do
        do k = 1, (n - 1)/2
            z(n - k + 1) = conjg(c(k + 1))
        end do
        call fft_c2c(z, 1)
        z = z/real(n, dp)
        tol = tol_fac*real(n, dp)*maxval(abs(x))
        do k = 1, n
            if (abs(z(k) - x(k)) > tol) then
                call print_fail
                print *, "n = ", n, ", sample ", k, ": got ", z(k), ", expected ", x(k)
                stop "fft_r2c/fft_c2c roundtrip deviates from input"
            end if
        end do
    end subroutine check_roundtrip

    subroutine fill_random(x)
        real(dp), intent(out) :: x(:)
        integer, allocatable :: seed(:)
        integer :: nseed

        call random_seed(size=nseed)
        allocate (seed(nseed))
        seed = 20260609
        call random_seed(put=seed)
        call random_number(x)
        x = 2.0d0*x - 1.0d0
    end subroutine fill_random

end program test_fft
