program ql_balance

use grid_mod
use baseparam_mod
use control_mod
use wave_code_data
use recstep_mod, only : nstack,tol,tim_stack,y_stack,timstep_arr
!DIAG:
use diag_mod, only : write_diag,iunit_diag,write_diag_b,iunit_diag_b
!END DIAG

use mpi

implicit none

integer :: ierror,np_num,irank

logical :: toomuch,opnd,dostep,scratch
logical :: flag_run_time_evolution !Added by Philipp Ulbl 12.05.2020
integer :: npoimin,ipoi,i,nstep,nstepmax,Nstorage,npoi,k,ieq,l
integer :: nmult,istage,itrans,ntrans,iunit_redo,ioddeven
double precision :: evoltime, timescale,timstep,eps,tmax,timstep_rec
double precision :: tmax_factor,antenna_factor
double precision :: antenna_factor_max !Added by Philipp Ulbl 12.05.2020
double precision :: relchg,relchgmax,facdecr,timstepmax
double precision :: err_tot,err_loc,err_minfac,tol_redfac
double precision :: tol_min,epsnoise,w,timstep_red,tol_max
double precision :: urelax,time,factolmax,factolred,timstep_min
double precision :: timscal_dql,rate_dql,timscal_dqli,rate_dqli
double precision :: stop_time_step !Added by Philipp Ulbl 13.05.2020
DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: yprev
DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: dqle11_prev
DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: dqle12_prev
DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: dqle21_prev
DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: dqle22_prev
DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: dqli11_prev
DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: dqli12_prev
DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: dqli21_prev
DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: dqli22_prev
DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: ych_one,ych_tot
double precision, dimension(:),   allocatable :: timscal
double precision, dimension(:),   allocatable :: dummy
double precision, dimension(:,:), allocatable :: params_beg,params_begbeg
double precision, dimension(:,:), allocatable :: params_num,params_denom
integer, dimension(1) :: ind_dqle,ind_dqli
  
!needed for interpolation of br abs and stopping criterion
!Added by Philipp Ulbl 04.06.2020
DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: br_abs
integer ::  ibrabsres, ibeg, iend, nlagr, nder
double precision, dimension(:,:),   allocatable :: coef
  
call MPI_Init(ierror);
call MPI_Comm_size(MPI_COMM_WORLD,np_num,ierror);
call MPI_Comm_rank(MPI_COMM_WORLD,irank,ierror);

if (irank .eq. 0 ) then
    print *,' '
    print *, '******************************'
    print *, 'number of processes:', np_num
    print *, '              irank:', irank
    print *, '******************************'
endif

open(22,file='balance.in');
read(22,*) flre_path;
read(22,*) vac_path;
read(22,*) btor;
read(22,*) rtor;
read(22,*) rmin;
read(22,*) rmax;
read(22,*) rsepar;
read(22,*) npoimin;
read(22,*) gg_factor;
read(22,*) gg_width;
read(22,*) gg_r_res;
read(22,*) Nstorage;
read(22,*) tmax_factor;
read(22,*) antenna_factor;
read(22,*) iboutype;
read(22,*) iwrite;
read(22,*) eps;
read(22,*) dperp;
read(22,*) icoll;
read(22,*) Z_i;
read(22,*) am;
read(22,*) rb_cut_in
read(22,*) re_cut_in
read(22,*) rb_cut_out
read(22,*) re_cut_out
read(22,*) write_formfactors
read(22,*) flag_run_time_evolution
read(22,*) stop_time_step
close(22);

! reread data form KiLCA input files:
open(22,file=trim(flre_path)//'background.in');
read(22,*);
read(22,*) rtor;
read(22,*);
read(22,*) btor;
close(22);

relchgmax=0.1d0
facdecr=1.d1
urelax=0.5d0 !0.5d0  !0.9d0
nmult=1 !10
nstack=2
tol_max=3.d-2 !3.d-4 !3.d-3 !3.d-2
!err_minfac=0.1d0
err_minfac= 1.d-2 !2.0d0/sqrt(dfloat(nmult))
tol_redfac=0.5d0
tol_min=3.d-5 
epsnoise=1.d-8
factolmax=3.d0
factolred=0.5d0
ntrans=10
!
!mwind=100
mwind=10
!

if (irank .eq. 0 ) then
    print *, ''
    print *, 'balance code V7'
    print *, '====================================================================='
    print *, 'Run time evolution: ', flag_run_time_evolution
    print *, ''
    print *, 'Parameters from input file:'
    print *, 'flre path: ', trim(flre_path)
    print *, 'vac path: ', trim(vac_path)
    print *, 'B_tor = ', btor
    print *, 'R_tor = ', rtor
    print *, 'r_min = ', rmin
    print *, 'r_max = ', rmax
    print *, 'npoimin = ', npoimin
    print *, 'gg_factor = ', gg_factor
    print *, 'gg_width = ', gg_width
    print *, 'gg_r_res = ', gg_r_res
    print *, 'Nstorage = ', Nstorage
    print *, 'tmax_factor = ', tmax_factor
    print *, 'antenna_factor = ', antenna_factor
    print *, 'iboutype = ', iboutype
    print *, 'iwrite = ', iwrite
    print *, 'eps = ', eps
    print *, 'dperp = ', dperp
    print *, 'icoll = ', icoll
    print *, 'Z_i = ', Z_i
    print *, 'am = ', am
    print *, 'stop_time_step = ', stop_time_step
    print *, ''
endif

timescale = (rmax-rmin)**2/dperp
tmax = timescale*tmax_factor
timstep = tmax/Nstorage
if (irank .eq. 0 ) then
    print *, 'timstep = ', timstep
endif
timstepmax=tmax

call gengrid(npoimin)
print *, 'irank = , npoib = ', irank, npoib

if(iboutype.eq.1) then
    npoi=npoic-1
else
    npoi=npoic
endif

allocate(yprev(neqset))
allocate(dqle11_prev(npoib))
allocate(dqle12_prev(npoib))
allocate(dqle21_prev(npoib))
allocate(dqle22_prev(npoib))
allocate(dqli11_prev(npoib))
allocate(dqli12_prev(npoib))
allocate(dqli21_prev(npoib))
allocate(dqli22_prev(npoib))
allocate(timstep_arr(neqset),tim_stack(neqset))

call initialize_wave_code_interface(npoib, rb);

!initial background profiles:
do ipoi=1,npoic
    !safety factor:
    qsaf(ipoi) = 0.5*(q(ipoi)+q(ipoi+1))
    !electron density :
    params(1,ipoi) = 0.5*(n(ipoi)+n(ipoi+1))
    !toroidal rotation frequency :
    params(2,ipoi) = 0.5*(Vz(ipoi)+Vz(ipoi+1))/rtor
    !electron temeperature :
    params(3,ipoi) = 0.5*(Te(ipoi)+Te(ipoi+1))*ev
    !ion temeperature :
    params(4,ipoi) = 0.5*(Ti(ipoi)+Ti(ipoi+1))*ev
end do

if (irank .eq. 0 ) then
    open(123,form='unformatted',file='init_params.dat')
    write(123) params
    close(123)
end if

call geomparprof

!
irf = 2
call get_dql
if(flag_run_time_evolution) then
    !For time evolution mode use antenna_factor as maximum
    !and start with a very small value and ramp this up
    !Added by Philipp Ulbl 12.05.2020
    antenna_factor_max = antenna_factor
    antenna_factor = 1.d-4
end if
dqle11=dqle11*antenna_factor
dqle12=dqle12*antenna_factor
dqle21=dqle21*antenna_factor
dqle22=dqle22*antenna_factor
dqli11=dqli11*antenna_factor
dqli12=dqli12*antenna_factor
dqli21=dqli21*antenna_factor
dqli22=dqli22*antenna_factor
irf = 1
!
call genstartsource

if (irank .eq. 0 ) then
    !for debugging:
    do ipoi=1,npoic
        !params(:,ipoi)=params(:,npoic)+1.0d0*(params(:,ipoi)-params(:,npoic))
        write (1000,*) rc(ipoi), params(1:2,ipoi)   &
                , params(3,ipoi)/ev                 &
                , params(4,ipoi)/ev                 &
                , 0.5d0*(Ercov(ipoi)+Ercov(ipoi+1)) &
                , 0.5d0*(sqg_bthet_overc(ipoi)+sqg_bthet_overc(ipoi+1))
    end do
close(1000)
end if

!stop !Stop Schalter verhindert Ausfuehrung des DGL-Loesers 

allocate(timscal(npoi),dummy(npoic))
allocate(params_beg(nbaleqs,npoic),params_num(nbaleqs,npoic))
allocate(params_denom(nbaleqs,npoic))
allocate(params_begbeg(nbaleqs,npoic))
time=0.d0
itrans=0
tol=tol_max
istage=1
!DIAG:
write_diag=.false.
write_diag_b=.false.
!END DIAG
!
inquire(file='restart.dat',exist=opnd)
!
if(opnd) then
    !
    if (irank .eq. 0 ) then
        print *,'restart'
    endif
    !
    open(201,file='final.restart')
    do ipoi=1,npoi
        read(201,*) timstep,params(:,ipoi)
        params(3:4,ipoi)=params(3:4,ipoi)*ev
        do ieq=1,nbaleqs
            k=nbaleqs*(ipoi-1)+ieq
            y(k)=params(ieq,ipoi)
        enddo
    enddo
    close(201)
    open(201,file='restart.dat')
    read(201,*) timstep
    close(201)
    scratch=.false.
    !
else
    !
    if (irank .eq. 0 ) then
        print *,'start from scratch'
    endif
    timstep=timstep*tol
    timstep=1.0d-12
    scratch=.true.
    !
endif
!
timstep_arr=0.d0
!
!  call evolvestep(timstep, eps)
!
timstep_arr=timstep
tim_stack=timstep_arr
!
!
print *,'start balance, irank = ', irank
iunit_diag=5000
write_diag=.true.
!
!pause
call get_dql
dqle11=dqle11*antenna_factor
dqle12=dqle12*antenna_factor
dqle21=dqle21*antenna_factor
dqle22=dqle22*antenna_factor
dqli11=dqli11*antenna_factor
dqli12=dqli12*antenna_factor
dqli21=dqli21*antenna_factor
dqli22=dqli22*antenna_factor
!pause
!!call MPI_finalize(ierror);
!!stop
!
allocate(br_abs(Nstorage))
dqle11_prev=dqle11
dqle12_prev=dqle12
dqle21_prev=dqle21
dqle22_prev=dqle22
dqli11_prev=dqli11
dqli12_prev=dqli12
dqli21_prev=dqli21
dqli22_prev=dqli22
params_begbeg=params
if (irank .eq. 0 ) then
    print *,'dql ready'
endif
!
!call equipotentials
!if (irank .eq. 0 ) then
!    print *,'equipotentials ready'
!endif

if(.not. flag_run_time_evolution) then
    !Stop if mode is not time evolution
    !Added by Philipp Ulbl 12.05.2020

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    print *, 'stop: linear code only'
    call MPI_finalize(ierror);
    stop  !! <<----- Stop for linear code usage
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
endif

!init variables for interpolation of Br abs res
!Added by Philipp Ulbl 04.06.2020
nlagr = 4;
nder = 0;
allocate(coef(0:nder,nlagr))
        
iunit_diag=137
iunit_diag_b=8138
!write_diag_b=.true.
ioddeven=1
open(4321,file='timstep_evol.dat')
open(777,file='br_abs_res.dat')
close(777)
!
do i=1,Nstorage
    !
    do ipoi=1,npoi
        do ieq=1,nbaleqs
            k=nbaleqs*(ipoi-1)+ieq
            yprev(k)=params(ieq,ipoi)
        enddo
    enddo
    !
    !
    dostep=.true.
        !
    do while(dostep)

        if (irank .eq. 0 ) then
            !
            !DIAG:
            iunit_diag=5000+i
            if(write_diag_b) then
                if(ioddeven/2*2.eq.ioddeven) then
                    open(iunit_diag_b,file='params_b_redostep.even')
                else
                    open(iunit_diag_b,file='params_b_redostep.odd')
                endif
            endif
            !END DIAG
            !
        end if

        call get_dql
        
        !write fort.7000 (debug)
        !the difference between fort5000 and fort7000 turned out to be some compiling mistake
        !open(iunit_diag+2000)
        !do ipoi = 1, npoib
        !    write(iunit_diag+2000,*) r(ipoi),dqle11(ipoi),dqle12(ipoi) &
        !        ,dqle22(ipoi),dqli11(ipoi)        &
        !        ,dqli12(ipoi),dqli22(ipoi)        &
        !        ,abs(Br(ipoi))                    &
        !        ,abs(Br(ipoi)-c*kp(ipoi)*Es(ipoi)/om_E(ipoi)) &
        !        ,abs(Br(ipoi)-c*ks(ipoi)*Ep(ipoi)/om_E(ipoi)) &
        !        ,abs(Jpe(ipoi)),abs(Jpi(ipoi)) &
        !        ,abs(Jpe(ipoi)+Jpi(ipoi)) 
        !enddo
        !close(iunit_diag+2000)
        
        !Stop if timestep becomes too small
        !Added by Philipp Ulbl 13.05.2020
        if (timstep .lt. stop_time_step .and. time .gt. 1.0d-3) then
            print *, 'stop: timestep smaller than stop limit'
            call MPI_finalize(ierror);
            stop
        endif
        
        !calculate Br abs at the resonant surface for stopping criterion
        !Added by Philipp Ulbl 04.06.2020

        !binsearch
        call binsrc(rb,1,npoib,r_resonant,ibrabsres)
        !
        ibeg=max(1,ibrabsres-nlagr/2)
        iend=ibeg+nlagr-1
        if(iend.gt.npoib) then
            iend=npoib
            ibeg=iend-nlagr+1
        endif

        !lagrange interpolation with order 4 only for function (0)
        call plag_coeff(nlagr,nder,r_resonant,rb(ibeg:iend),coef)
        br_abs(i)=sum(coef(0,:)*abs(Br(ibeg:iend)))*sqrt(antenna_factor)
        
        !output on console and save to file
        print *, 'Br abs res = ', br_abs(i)
        open(777,file='br_abs_res.dat',position='append')
        write(777,*) i,time,antenna_factor,br_abs(i)
        close(777)
        
        !Stop if gradient of Abs(Br) at the resonance becomes negative
        if (i .gt. 5.d0 .and. br_abs(i)-br_abs(i-1) .lt. 0.d0 &
                        .and. br_abs(i-1)-br_abs(i-2) .lt. 0.d0) then
            print *, 'stop: gradient of Br_Abs_Res negative'
            !call MPI_finalize(ierror);
            !stop
        endif
        
        !Ramp up antenna_factor: linear in Icoil or quadratic in D
        !Added by Philipp Ulbl 12.05.2020
        if (antenna_factor .lt. (5.d0*antenna_factor_max)) then !5x for testing purposes
            antenna_factor = time**2 + 1.d-4
            !This can be activated for runs without QL evolution to check steady state behaviour
            !antenna_factor = 1.d-4
            !if(i .gt. 200) then
            !    stop
            !endif
            print *,'antenna_factor = ',antenna_factor
        else
            print *, 'stop: reached antenna_factor_max'
            call MPI_finalize(ierror);
            stop
        endif

        !Old ramp up by Martin: ramp up to 1s then ramp down to 2s
        !to check hysteresis behaviour. afterwards small increases
        !if (time.lt.1.d0) then
        !    antenna_factor = time**2
        !elseif (time.gt.1.d0.AND.time.lt.2.d0 ) then
        !    antenna_factor = (2.d0 - time)**2
        !endif
        !antenna_factor = antenna_factor + 1.d-4

        dqle11=dqle11*antenna_factor
        dqle12=dqle12*antenna_factor
        dqle21=dqle21*antenna_factor
        dqle22=dqle22*antenna_factor
        dqli11=dqli11*antenna_factor
        dqli12=dqli12*antenna_factor
        dqli21=dqli21*antenna_factor
        dqli22=dqli22*antenna_factor

        !
        !DIAG:
        if(write_diag_b) close(iunit_diag_b)
        !END DIAG
        timscal_dql=maxval(abs(dqle11_prev-dqle11))/maxval(dqle11_prev+dqle11)
        ind_dqle=maxloc(abs(dqle11_prev-dqle11))
        timscal_dqli=maxval(abs(dqli11_prev-dqli11))/maxval(dqli11_prev+dqli11)
        ind_dqli=maxloc(abs(dqli11_prev-dqli11))
        rate_dql=timscal_dql/timstep
        rate_dqli=timscal_dqli/timstep
        if (irank .eq. 0 ) then
            print *,'timscal_dqle = ',sngl(timscal_dql) &
                    ,'timscal_dqli = ',sngl(timscal_dqli)
            print *,'maximum dqle at r = ',rc(ind_dqle(1)) &
                    ,'maximum dqli at r = ',rc(ind_dqli(1))
            do ipoi=1,npoib
                write(9999,*) rb(ipoi),abs(dqle11_prev(ipoi)-dqle11(ipoi)),abs(dqli11_prev(ipoi)-dqli11(ipoi))
            enddo
            close(9999)
        endif
        !    timscal_dql=timscal_dql+timscal_dqli
        !    if(timscal_dql.lt.tol*factolmax) then
        if(.true.) then
            dostep=.false.
            dqle11_prev=dqle11
            dqle12_prev=dqle12
            dqle21_prev=dqle21
            dqle22_prev=dqle22
            dqli11_prev=dqli11
            dqli12_prev=dqli12
            dqli21_prev=dqli21
            dqli22_prev=dqli22
            params_begbeg=params
            ioddeven=ioddeven+1
        else
            if (irank .eq. 0 ) then
                print *,'redo step with old DQL'
            endif
            dqle11=dqle11_prev
            dqle12=dqle12_prev
            dqle21=dqle21_prev
            dqle22=dqle22_prev
            dqli11=dqli11_prev
            dqli12=dqli12_prev
            dqli21=dqli21_prev
            dqli22=dqli22_prev
            iunit_redo=137
            if (irank .eq. 0 ) then
                open(iunit_redo,file='params_redostep.after')
                do ipoi=1,npoic
                    write (iunit_redo,*) rc(ipoi),params(1:2,ipoi)    &
                            , params(3,ipoi)/ev                  &
                            , params(4,ipoi)/ev                  &
                            , 0.5d0*(Ercov(ipoi)+Ercov(ipoi+1))
                end do
                close(iunit_redo)
            end if
            params=params_begbeg
            if (irank .eq. 0 ) then
                open(iunit_redo,file='params_redostep.before')
                do ipoi=1,npoic
                    write (iunit_redo,*) rc(ipoi),params(1:2,ipoi)    &
                            , params(3,ipoi)/ev                  &
                            , params(4,ipoi)/ev                  &
                            , 0.5d0*(Ercov(ipoi)+Ercov(ipoi+1))
                end do
                close(iunit_redo)
            end if
            timstep=timstep/factolmax
            timstep_arr=timstep
        endif
        !
        do
        !
            params_beg=params
        !
            call evolvestep(timstep, eps)
        !
            !limits ion and electron temperatures from below (by 10 eV in this example).
            !Added by Philipp Ulbl on 09.06.2020
            do ipoi=1,npoic
                params(3,ipoi)=max(params(3,ipoi),50d0*ev)
                params(4,ipoi)=max(params(4,ipoi),50d0*ev)
            enddo 
        !
            params_num=(params-params_beg)**2
            params_denom=params**2+params_beg**2
        !
            do ieq=1,nbaleqs
        !
                call smooth_array_gauss(npoi,mwind,params_num(ieq,:),dummy)
        !
                params_num(ieq,:)=dummy
        !
                call smooth_array_gauss(npoi,mwind,params_denom(ieq,:),dummy)
        !
                params_denom(ieq,:)=dummy
            enddo 
        !
            do ipoi=1,npoi
                if(rc(ipoi).lt.0.95d0*rc(npoi)) then
                    timscal(ipoi)=sum(sqrt(params_num(3:4,ipoi)/params_denom(3:4,ipoi)))
                else
                    timscal(ipoi)=1d-30 !0.d0
                endif
            enddo
        !
            if (irank .eq. 0 ) then
                print *,'maxval(timscal) = ',maxval(timscal)
            endif
            if(maxval(timscal).lt.tol*factolmax) exit
        !
            timstep_arr=timstep_arr*factolred
            params=params_beg
            if (irank .eq. 0 ) then 
                print *,'redo step'
            endif
        enddo
    !
        timscal=timscal+timscal_dql
    !
        do ipoi=1,npoi
            do ieq=1,nbaleqs
                k=nbaleqs*(ipoi-1)+ieq
                !timstep_arr(k)=timstep_arr(k)/timscal(ipoi)*tol
                timstep_arr(k)=timstep_arr(k)/max(timscal(ipoi),epsilon(1.d0))*tol
            enddo
        enddo
    !
        timstep_arr=timstep_arr*timescale/(timstep_arr+timescale)
        if(scratch) then
            scratch=.false.
            tim_stack=timstep_arr
        endif
        timstep_arr=2.d0*timstep_arr*tim_stack/(timstep_arr+tim_stack)
        timstep=minval(timstep_arr)
        if(time .gt. 1.2 .and. timstep .gt. 1.d-4) then
            timstep = 1.d-4
        endif
        open(5432,file='timstep_min.inp')
        read (5432,*) timstep_min
        close(5432)
        timstep=max(timstep,timstep_min)
        !
        timstep_arr=timstep
        !
        tim_stack=timstep_arr
        !
        if (irank .eq. 0 ) then
            print *,'timstep',real(timstep),'   timescale',real(timescale), &
                'tolerance',real(tol)
        endif
    !
        if (irank .eq. 0 ) then
            write(4321,*) i,timstep,timscal_dql,timscal(1),rate_dql,time
            close(4321)
            open(4321,file='timstep_evol.dat',position='append')
        end if
    !
    enddo
    !
    do ipoi=1,npoi
        do ieq=1,nbaleqs
            k=nbaleqs*(ipoi-1)+ieq
            params(ieq,ipoi)=yprev(k)*urelax+params(ieq,ipoi)*(1.d0-urelax)
        enddo
    enddo
    timstep_arr=0.d0
    call evolvestep(timstep, eps)
    timstep_arr=timstep
    !
    time=time+timstep
    !   
    !
    if (irank .eq. 0 ) then
        print *, 'i = ', int2(i), 'time = ', real(time)
        print *,' '
    endif
    !
    !for debugging:
    if (irank .eq. 0 ) then
        do ipoi=1,npoic
            write (1000+i,*) rc(ipoi),params(1:2,ipoi)    &
                    , params(3,ipoi)/ev                  &
                    , params(4,ipoi)/ev                  &
                    , 0.5d0*(Ercov(ipoi)+Ercov(ipoi+1))  &
                    , 0.5d0*(sqg_bthet_overc(ipoi)+sqg_bthet_overc(ipoi+1))
        end do
        close(1000+i)
    end if
!
end do

call deallocate_wave_code_data ();
deallocate(yprev)
deallocate(tim_stack)
deallocate(timscal,params_beg,params_num,params_denom,dummy)

print *, 'Programm is finalized';
call MPI_finalize(ierror);

end
