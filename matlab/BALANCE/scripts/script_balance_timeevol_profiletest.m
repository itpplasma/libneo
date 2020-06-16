%##########################################################################
% script_experimental.m
%##########################################################################
% description:
%--------------------------------------------------------------------------
% this script is for experimental tries with the balance code.
%##########################################################################

%author:   Philipp Ulbl
%created:  05.05.2020

libKiLCA = '~/KiLCA_interface/';
libBalance = '~/BALANCE/balance';

addpath(genpath(libKiLCA))
addpath(genpath(libBalance))

mpath = pwd();

addpath('~/BALANCE/balance');

studyname = 'TimeEvol_ProfTest';
system(['mkdir -p ~/Balance_Results/', studyname, '/']);

shot = 33133;
time = 3000;

time_evol = true;

m = 5;
n = 2 .* ones(size(m));

%input
gfile  = ['/temp/ulbl_p/AUG/SHOTS/',num2str(shot),'/g',num2str(shot),'.',num2str(time),'_EQH'];
filehead = ['/temp/ulbl_p/AUG/SHOTS/',num2str(shot),'/',num2str(shot),'.',num2str(time)];
cfile  = [filehead,'_coil.dat'];

dapath = '/temp/ulbl_p/DA/ASTRA/';
%dapath = '';

neprof = [filehead,'_ne_PED_ULBLP_rho_pol.dat'];
Teprof = [filehead,'_Te_PED_ULBLP_rho_pol.dat'];
Tiprof = [filehead,'_Ti_PED_ULBLP_rho_pol.dat'];
vtprof = [filehead,'_vt_PED_ULBLP_rho_pol.dat'];

%location of field.dat
%pfile = ['/temp/ulbl_p/MESH3D/',num2str(shot),'/field.dat']; %will be calculated if not present
fluxdatapath = ['/temp/ulbl_p/FLUXDATA/',num2str(shot),'/',num2str(time),'/']; %will be calculated if not present

gpecpath = ['/temp/ulbl_p/GPEC/TimeEvol/', num2str(shot), '_', num2str(time),'/'];
%copy = '/temp/ulbl_p/BALANCE_2020/TimeEvol_m6_experimental/33133_3000/profiles/';
%copy='';

dname = '/temp/ulbl_p/PROFILES/33133_3000_TimeEvolTest/';
profs = dir(dname);
profs = profs(3:end); %remove . and ..

for k = 4:6%2:numel(profs)

    copy = [dname, profs(k).name, '/'];
    sname = [studyname, '_', profs(k).name];
    runpath = ['/temp/ulbl_p/BALANCE_2020/', studyname, '/', profs(k).name, '/', num2str(shot), '_', num2str(time),'/'];

    disp(sname);
    
    %##########################################################################
    % RUN BALANCE WITH MATLAB CLASS
    %##########################################################################

    %BALANCE CODE
    bal = Balance(runpath, shot, time, sname);
    bal.setModes(m, n);
    bal.setCoil(cfile);
    bal.setEqui(gfile, fluxdatapath);
    bal.setTMHDCode('GPEC', gpecpath);
    bal.setProfiles(neprof, Teprof, Tiprof, vtprof, copy);
    bal.setKiLCA();
    bal.setDaEstimation(dapath);
    opt = balanceoptions(bal.kil_flre.pathofrun, bal.kil_vacuum.pathofrun);
    opt.stop_time_step=1e-8;
    opt.flag_run_time_evolution = time_evol;
    bal.setOptions(opt);
    bal.write();
    bal.run();

    %##########################################################################
    % EXPORT 2 HDF5 FORMAT
    %##########################################################################

    system(['mkdir -p ~/Balance_Results/', studyname, '/']);
    bal.export2HDF5(['~/Balance_Results/', studyname, '/'], sname);

end









