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

studyname = 'TimeEvolution_TestKiLCABack';
system(['mkdir -p ~/Balance_Results/', studyname, '/']);

shot = 33133;
time = 3000;

time_evol = true;

m = 5;
n = 2 .* ones(size(m));

runpath = ['/temp/ulbl_p/BALANCE_2020/', studyname, '/', num2str(shot), '_', num2str(time),'/'];

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
copy='/temp/ulbl_p/BALANCE_2020/TimeEvolution/33133_3000/profiles/';

%REF FROM MARTIN
% gfile  = '/proj/plasma/RMP/DATA2017/33133/3.0s/g33133.3000_ed4';
% cfile  = '/temp/ulbl_p/AUG/SHOTS/33120/33120.5500_coil.dat';
% 
% neprof = '/proj/plasma/RMP/DATA2017/33133/3.0s/neprof_aug_33133_3.asc';
% Teprof = '/proj/plasma/RMP/DATA2017/33133/3.0s/Teprof_aug_33133_3.asc';
% Tiprof = '/proj/plasma/RMP/DATA2017/33133/3.0s/Tiprof_aug_33133_3.asc';
% vtprof = '/proj/plasma/RMP/DATA2017/33133/3.0s/vtprof_aug_33133_3.asc';
% 
% dapath = '/temp/ulbl_p/DA/ASTRA/';
% %pfile = '/temp/kasilov/MESH3D_33133/field.dat';
% %fluxdatapath = '/proj/plasma/RMP/DATA2017/33133/3.0s/PREPROCESS/';
% %gpecpath = ['/temp/ulbl_p/GPEC/TimeEvol_ReprMartin33133/', num2str(shot), '_', num2str(time),'/'];
% copy = '/temp/ulbl_p/BALANCE_2020/TimeEvol_ReprMartin33133/33133_3000/profiles/';
% 
% fluxdatapath = ['/temp/ulbl_p/FLUXDATA/',num2str(shot),'/',num2str(time),'/']; %will be calculated if not present
% gpecpath = ['/temp/ulbl_p/GPEC/TimeEvol/', num2str(shot), '_', num2str(time),'/'];

%##########################################################################
% RUN BALANCE WITH MATLAB CLASS
%##########################################################################

%BALANCE CODE
bal = Balance(runpath, shot, time, studyname);
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
bal.kil_flre.background.flag_recalc = -1;
bal.kil_vacuum.background.flag_recalc = -1;
bal.kil_flre.write();
bal.kil_vacuum.write();
bal.run();

%##########################################################################
% EXPORT 2 HDF5 FORMAT
%##########################################################################

system(['mkdir -p ~/Balance_Results/', studyname, '/']);
bal.export2HDF5(['~/Balance_Results/', studyname, '/'], studyname);











