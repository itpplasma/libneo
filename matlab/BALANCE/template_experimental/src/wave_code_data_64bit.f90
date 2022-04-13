!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module wave_code_data

implicit none;

integer, parameter :: pp = 8;   !pp = 4 for 32 bit and pp = 8 for 64 bit

integer(pp), allocatable, dimension(:) :: flre_cd_ptr; !pointer to kilca's core data object for each mode
integer(pp), allocatable, dimension(:) :: vac_cd_ptr;  !pointer to kilca's core data object for each mode

character(1024) :: flre_path;
character(1024) :: vac_path;

integer :: flre_call_ind = 0;
integer :: vac_call_ind = 0;

integer :: dim_mn;
integer, allocatable, dimension(:) :: m_vals, n_vals;

!radial grid:
integer :: dim_r;                            !radial grid dimension
real(8), allocatable, dimension(:) :: r;     !radial grid array

!wave fields (Gauss system units) for a mode:
complex(8), allocatable, dimension(:) :: Er; !E field in rsp coordinates
complex(8), allocatable, dimension(:) :: Es; !E field in rsp coordinates
complex(8), allocatable, dimension(:) :: Ep; !E field in rsp coordinates
complex(8), allocatable, dimension(:) :: Et; !Etheta field in cyl coordinates
complex(8), allocatable, dimension(:) :: Ez; !Ez field in cyl coordinates

complex(8), allocatable, dimension(:) :: Br; !B field in rsp coordinates
complex(8), allocatable, dimension(:) :: Bs; !B field in rsp coordinates
complex(8), allocatable, dimension(:) :: Bp; !B field in rsp coordinates
complex(8), allocatable, dimension(:) :: Bt; !Btheta field in cyl coordinates
complex(8), allocatable, dimension(:) :: Bz; !Bz field in cyl coordinates

!current densities (Gauss system units) for a mode:
complex(8), allocatable, dimension(:) :: Jri; !Jr for ions
complex(8), allocatable, dimension(:) :: Jsi; !Js for ions
complex(8), allocatable, dimension(:) :: Jpi; !Jp for ions
complex(8), allocatable, dimension(:) :: Jre; !Jr for electrons
complex(8), allocatable, dimension(:) :: Jse; !Js for electrons
complex(8), allocatable, dimension(:) :: Jpe; !Jp for electrons

!background profiles:
real(8), allocatable, dimension(:) :: q;     !safety factor
real(8), allocatable, dimension(:) :: n;     !1/cm^3
real(8), allocatable, dimension(:) :: Ti;    !eV
real(8), allocatable, dimension(:) :: Te;    !eV
real(8), allocatable, dimension(:) :: Vth;   !cm/c
real(8), allocatable, dimension(:) :: Vz;    !cm/c
real(8), allocatable, dimension(:) :: dPhi0; !electric potential, Gauss units

real(8), allocatable, dimension(:) :: nui;   !ions collision frequency
real(8), allocatable, dimension(:) :: nue;   !electrons collision frequency

real(8), allocatable, dimension(:) :: B0t;
real(8), allocatable, dimension(:) :: B0z;
real(8), allocatable, dimension(:) :: B0;

!misc data for a mode:
real(8), allocatable, dimension(:) :: kp;
real(8), allocatable, dimension(:) :: ks;

real(8), allocatable, dimension(:) :: om_E;

!for a spectrum:
real(8), allocatable, dimension(:) :: diss_pow_dens;

! initial background profiles to be loaded from files:
real(8), allocatable, dimension(:) :: rq,   iq;     !safety factor
real(8), allocatable, dimension(:) :: rn,   in;     !1/cm^3
real(8), allocatable, dimension(:) :: rTi,  iTi;    !eV
real(8), allocatable, dimension(:) :: rTe,  iTe;    !eV
real(8), allocatable, dimension(:) :: rVth, iVth;   !cm/c
real(8), allocatable, dimension(:) :: rVz,  iVz;    !cm/c
real(8), allocatable, dimension(:) :: rep,  idPhi0; !electric potential, Gauss units

character(1024) :: path2profs = './profiles/';      ! path to background profiles

end module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_wave_code_data(imin, imax)

!runs wave code to compute (E,B) - fields.

use wave_code_data;

implicit none;

integer, intent (in) :: imin, imax;

integer :: k;

do k = imin,imax

    call clear_wave_code_data(flre_cd_ptr(k));

    call calc_wave_code_data_for_mode(flre_cd_ptr(k), flre_path, len(trim(flre_path)), m_vals(k), n_vals(k));

end do

flre_call_ind = flre_call_ind + 1;

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine initialize_wave_code_interface(nrad, r_grid)

use wave_code_data;

implicit none;

integer, intent (in) :: nrad;
real(8), dimension(nrad), intent (in) :: r_grid;

integer :: k;

call read_antenna_modes(flre_path); ! read antenna modes from modes.in file and allocate arrays for mode's numbers

call allocate_wave_code_data(nrad, r_grid);

call read_background_profiles(path2profs);

call interp_background_profiles(); ! interpolate the initial profiles to the balance grid

! vacuum fields for the whole spectrum:
do k = 1,dim_mn

    call clear_wave_code_data(vac_cd_ptr(k));

    call calc_wave_code_data_for_mode(vac_cd_ptr(k), vac_path, len(trim(vac_path)), m_vals(k), n_vals(k));

end do

call get_background_magnetic_fields_from_wave_code(vac_cd_ptr(1), dim_r, r, B0t, B0z, B0);

call get_collision_frequences_from_wave_code(vac_cd_ptr(1), dim_r, r, nui, nue);

vac_call_ind = vac_call_ind + 1;

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine allocate_wave_code_data (d_grid, r_grid)

use wave_code_data;

implicit none;

integer, intent (in) :: d_grid;
real(8), dimension(d_grid), intent (in) :: r_grid;

dim_r = d_grid;

allocate (r(dim_r));

allocate (Er(dim_r), Es(dim_r), Ep(dim_r), Et(dim_r), Ez(dim_r));

allocate (Br(dim_r), Bs(dim_r), Bp(dim_r), Bt(dim_r), Bz(dim_r));

allocate (Jri(dim_r), Jsi(dim_r), Jpi(dim_r), Jre(dim_r), Jse(dim_r), Jpe(dim_r));

allocate (q(dim_r), n(dim_r), Ti(dim_r), Te(dim_r), Vth(dim_r), Vz(dim_r), dPhi0(dim_r));

allocate (kp(dim_r), ks(dim_r));

allocate (om_E(dim_r));

allocate (B0t(dim_r), B0z(dim_r), B0(dim_r));

allocate (nui(dim_r), nue(dim_r));

allocate (diss_pow_dens(dim_r));

r = r_grid;

allocate(flre_cd_ptr(dim_mn), vac_cd_ptr(dim_mn));

flre_cd_ptr = 0;
vac_cd_ptr = 0;

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine deallocate_wave_code_data ()

use wave_code_data;

implicit none;

integer :: k;

!removes wave_code objects:
do k = 1,dim_mn

    call clear_wave_code_data (vac_cd_ptr(k));
    call clear_wave_code_data (flre_cd_ptr(k));

end do

deallocate (r);

deallocate (Er, Es, Ep, Et, Ez);

deallocate (Br, Bs, Bp, Bt, Bz);

deallocate (Jri, Jsi, Jpi, Jre, Jse, Jpe);

deallocate (q, n, Ti, Te, Vth, Vz, dPhi0);

deallocate (kp, ks);

deallocate (om_E);

deallocate (B0t, B0z, B0);

deallocate (nui, nue);

deallocate (diss_pow_dens);

deallocate (m_vals, n_vals);

deallocate(flre_cd_ptr, vac_cd_ptr);

deallocate(rq,   iq);
deallocate(rn,   in);
deallocate(rTi,  iTi);
deallocate(rTe,  iTe);
deallocate(rVth, iVth);
deallocate(rVz,  iVz);
deallocate(rep,  idPhi0);

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine save_wave_code_data ()

use wave_code_data;

implicit none;

integer :: k;

character(64) :: string;
character(16) :: ind_str;

integer :: width;

if (flre_call_ind > 0) then
    width = log10(real(flre_call_ind)) + 1;
else
    width = 1;
end if

write(unit=ind_str, fmt='(I8)') width;

write(unit=ind_str, fmt=*) trim(ind_str);

write(unit=string, fmt='(A5,I'//ind_str//',A4)') 'back_', flre_call_ind, '.dat';
open (10, file=trim(string));

write(unit=string, fmt='(A5,I'//ind_str//',A4)') 'misc_', flre_call_ind, '.dat';
open (20, file=trim(string));

write(unit=string, fmt='(A7,I'//ind_str//',A4)') 'fields_', flre_call_ind, '.dat';
open (30, file=trim(string));

!write(unit=string, fmt='(A12,I'//ind_str//',A4)') 'disspowdens_', flre_call_ind, '.dat';
!open (40, file=trim(string));

do k=1,dim_r

    write (10,*) r(k), q(k), n(k), Ti(k), Te(k), Vth(k), Vz(k), dPhi0(k);

    write (20,*) r(k), kp(k), ks(k), B0t(k), B0z(k), B0(k), om_E(k);

    write (30,*) r(k), real(Er(k)), aimag(Er(k)), real(Es(k)), aimag(Es(k)), &
                       real(Ep(k)), aimag(Ep(k)), real(Et(k)), aimag(Et(k)), &
                       real(Ez(k)), aimag(Ez(k)), real(Br(k)), aimag(Br(k)), &
                       real(Bs(k)), aimag(Bs(k)), real(Bp(k)), aimag(Bp(k)), &
                       real(Bt(k)), aimag(Bt(k)), real(Bz(k)), aimag(Bz(k));

!    write (40,*) r(k), diss_pow_dens(k);

end do

close (10);
close (20);
close (30);
!close (40);

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine update_background_files (path)

use wave_code_data;
use grid_mod, only: params_b, Ercov;
use baseparam_mod;

implicit none;

character(1024), intent(in) :: path;

integer :: k;

integer :: flag = 0; ! set it 1 if you want to store the background profiles to disk

do k = 1,dim_r !at cell boundaries

    n(k)  = params_b(1,k);

    Te(k) = params_b(3,k)/ev;

    Ti(k) = params_b(4,k)/ev;

    Vz(k) = params_b(2,k)*rtor;

    dPhi0(k) = - Ercov(k);

end do

if (flag > 0) then

    open (10, file=trim(path)//'q.dat');
    open (20, file=trim(path)//'n.dat');
    open (30, file=trim(path)//'Te.dat');
    open (40, file=trim(path)//'Ti.dat');
    open (50, file=trim(path)//'Vth.dat');
    open (60, file=trim(path)//'Vz.dat');
    open (70, file=trim(path)//'Er.dat');

    do k=1,dim_r !at cell boundaries
        write (10,*) r(k), q(k);
        write (20,*) r(k), n(k);
        write (30,*) r(k), Te(k);
        write (40,*) r(k), Ti(k);
        write (50,*) r(k), Vth(k);
        write (60,*) r(k), Vz(k);
        write (70,*) r(k), - dPhi0(k);
    end do

    close (10);
    close (20);
    close (30);
    close (40);
    close (50);
    close (60);
    close (70);

end if

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_background_dimension_from_balance (dim_p)

use wave_code_data;

implicit none;

integer, intent(out) :: dim_p;

dim_p = dim_r;

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_background_profiles_from_balance (dim_p, r_p, q_p, n_p, Ti_p, Te_p, Vth_p, Vz_p, Er_p)

use wave_code_data;

implicit none;

integer, intent(out) :: dim_p;
real(8), dimension(dim_r), intent(out) :: r_p;
real(8), dimension(dim_r), intent(out) :: q_p;
real(8), dimension(dim_r), intent(out) :: n_p;
real(8), dimension(dim_r), intent(out) :: Ti_p;
real(8), dimension(dim_r), intent(out) :: Te_p;
real(8), dimension(dim_r), intent(out) :: Vth_p;
real(8), dimension(dim_r), intent(out) :: Vz_p;
real(8), dimension(dim_r), intent(out) :: Er_p;

dim_p = dim_r;
r_p   = r;
q_p   = q;
n_p   = n;
Ti_p  = Ti;
Te_p  = Te;
Vth_p = Vth;
Vz_p  = Vz;
Er_p  = - dPhi0;

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine eval_diss_power_density (dim, r, type, spec, d_p_d)

use wave_code_data, only: dim_mn, m_vals, n_vals, vac_cd_ptr, flre_cd_ptr;
use baseparam_mod, only: rtor;

implicit none;

integer, intent (in) :: dim;
real(8), dimension(dim), intent (in) :: r;
integer, intent (in) :: type, spec;
real(8), dimension(dim), intent (out) :: d_p_d;

complex(8), dimension(dim) :: amn_psi, amn_theta, amn_theta_cyl;

complex(8), dimension(dim) :: F, Br; !for vacuum fields

real(8), dimension(dim) :: dpd_mn;

integer :: ind, k, ierr;

d_p_d = 0.0d0;

do ind = 1,dim_mn !over modes

    !amn
    do k = 1,dim
        call amn_of_r (m_vals(ind), n_vals(ind), r(k), amn_psi(k), amn_theta(k), ierr);

        !if (ierr /= 0) then
        !    print *, 'amn_of_r: warning: ierr =', ierr, m_vals(ind), n_vals(ind), r(k);
        !end if
    end do

    call get_wave_fields_from_wave_code (vac_cd_ptr, dim, r, m_vals(ind), n_vals(ind), F, F, F, F, F, Br, F, F, F, F);

    call get_diss_power_density_from_wave_code (flre_cd_ptr, dim, r, m_vals(ind), n_vals(ind), type, spec, dpd_mn);

    amn_theta_cyl = (r*rtor/n_vals(ind)) * Br;

    d_p_d = d_p_d + 2.0d0 * dpd_mn * (abs(amn_theta)**2 / abs(amn_theta_cyl)**2);

end do

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine set_flag_for_profiles_in_background_input_file(path, flag)

character(128), intent(in) :: path;
integer, intent(in) :: flag;

character(len=1024) :: buffer;

character(len=1024), dimension(64) :: lines;

integer :: ios, num, ind;

open (10, file=trim(path)//'background.in');

ios = 0;
num = 0;

do while (ios == 0)

    read(10, '(A)', iostat=ios) buffer;

    if (ios == 0) then

        num = num + 1;

        if (num > 64) then
            write(*,*) 'Bad file format of background.in file!';
            exit;
        end if;

        lines(num) = buffer;

    end if

end do

close(10);

!modify the line with background profiles flag:
if      (flag == 1) then
    lines(8)(1:2) = '1 ';
else if (flag == 2) then
    lines(8)(1:2) = '2 ';
else if (flag == -1) then
    lines(8)(1:2) = '-1';
else if (flag == -2) then
    lines(8)(1:2) = '-2';
else
    write(*,*) 'Bad value of background profiles flag!';
end if

open (20, file=trim(path)//'background.in');

do ind = 1,num

    write(20, '(A)', iostat=ios) trim(lines(ind));

end do

close(20);

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine find_file_length(filename, l)

! find how many lines are in the filename

implicit none;

character(1024) :: filename; ! attention: pass here a filename with exactly the same length
integer :: l;

character(1024) :: buffer;
integer :: ios;

open (10, file=trim(filename));

ios = 0;
l = 0;

do while (ios == 0)

    read(10, '(A)', iostat=ios) buffer;

    if (ios == 0) then

        l = l + 1;

    end if

end do

close(10);

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine load_profile(filename, l, r, q)

! load a profile of length l nodes to pre-allocated arrays r, q

implicit none;

character(1024)       :: filename; ! attention: pass here a filename with exactly the same length
integer               :: l;
real(8), dimension(l) :: r;
real(8), dimension(l) :: q;

integer :: i;

open (10, file=trim(filename));

do i = 1,l

    read(10, *) r(i), q(i);

end do

close(10);

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine read_antenna_modes(path)

use wave_code_data, only: dim_mn, m_vals, n_vals;

implicit none;

character(1024) :: path;

character(1024) :: buffer;
integer :: l, i1, i2, i3;

open(10,file=trim(path)//'antenna.in');
read(10,*) buffer;
read(10,*) buffer;
read(10,*) buffer;
read(10,*) buffer;
read(10,*) buffer;
read(10,*) dim_mn;
close(10);

allocate (m_vals(dim_mn), n_vals(dim_mn));

open (10, file=trim(path)//'modes.in');

do l = 1,dim_mn

    read(10, '(A)') buffer;

    i1 = index(buffer, '(');
    i2 = index(buffer, ',');
    i3 = index(buffer, ')');

    read(buffer(i1+1:i2-1), *) m_vals(l);
    read(buffer(i2+1:i3-1), *) n_vals(l);

end do

close(10);

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine read_background_profiles(path)

use wave_code_data, only: rq, iq, rn, in, rTi, iTi, rTe, iTe, rVth, iVth, rVz, iVz, rep, idPhi0;

implicit none;

character(1024) :: path;
character(1024) :: file;
integer :: l;

file = trim(path)//'q.dat';
call find_file_length(file, l);
allocate(rq(l), iq(l));
call load_profile(file, l, rq, iq);

file = trim(path)//'n.dat';
call find_file_length(file, l);
allocate(rn(l), in(l));
call load_profile(file, l, rn, in);

file = trim(path)//'Ti.dat';
call find_file_length(file, l);
allocate(rTi(l), iTi(l));
call load_profile(file, l, rTi, iTi);

file = trim(path)//'Te.dat';
call find_file_length(file, l);
allocate(rTe(l), iTe(l));
call load_profile(file, l, rTe, iTe);

file = trim(path)//'Vth.dat';
call find_file_length(file, l);
allocate(rVth(l), iVth(l));
call load_profile(file, l, rVth, iVth);

file = trim(path)//'Vz.dat';
call find_file_length(file, l);
allocate(rVz(l), iVz(l));
call load_profile(file, l, rVz, iVz);

file = trim(path)//'Er.dat';
call find_file_length(file, l);
allocate(rep(l), idPhi0(l));
call load_profile(file, l, rep, idPhi0);

idPhi0 = - idPhi0; ! Er was loaded from Er.dat

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine interp_background_profiles()

! interpolate the initial profiles to the balance grid

use wave_code_data;

implicit none;

call interp_profile(size(rq),   rq,   iq,     dim_r, r, q);
call interp_profile(size(rn),   rn,   in,     dim_r, r, n);
call interp_profile(size(rTi),  rTi,  iTi,    dim_r, r, Ti);
call interp_profile(size(rTe),  rTe,  iTe,    dim_r, r, Te);
call interp_profile(size(rVth), rVth, iVth,   dim_r, r, Vth);
call interp_profile(size(rVz),  rVz,  iVz,    dim_r, r, Vz);
call interp_profile(size(rep),  rep,  idPhi0, dim_r, r, dPhi0);

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine interp_profile(dim_old, r_old, q_old, dim_new, r_new, q_new)

! interpolate a profile to a new grid by moving polynom of degree deg

implicit none;

integer :: dim_old;
real(8), dimension(dim_old) :: r_old, q_old;

integer :: dim_new;
real(8), dimension(dim_new) :: r_new, q_new;

integer :: deg = 9, Dmin = 0, Dmax = 0, l, ind;

do l = 1,dim_new

    call eval_neville_polynom(dim_old, r_old, q_old, deg, r_new(l), Dmin, Dmax, ind, q_new(l));

end do

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
