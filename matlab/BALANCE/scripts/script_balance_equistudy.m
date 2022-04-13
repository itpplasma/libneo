%##########################################################################
% script_balance_equistudy.m
%##########################################################################
% description:
%--------------------------------------------------------------------------
% this script runs the balance code for the study of the effect of differen
% t equilibria.
%##########################################################################

%author:   Philipp Ulbl
%created:  08.05.2020

libKiLCA = '~/KiLCA_interface/';
libBalance = '~/BALANCE/balance';

addpath(genpath(libKiLCA))
addpath(genpath(libBalance))
addpath(genpath('~/GITHUB/libneo/matlab/GPEC_interface/'));

mpath = pwd();

FLAG_FORCE_GPEC = false;
FLAG_FORCE_PROFILES = false;
FLAG_FORCE_NEO2 = false;
FLAG_FORCE_TMHD = true;

studyname = 'EquiStudy';
system(['mkdir -p ~/Balance_Results/', studyname, '/']);

%Runs to make
equis = {'33353_EQH',...
         '33353_AUGDIDEed0', ...
         '33353_MICDUEQBed0', ...
         '33353_MICDUEQBed1', ...
         '33353_MICDUEQBed2', ...
         %'33353_MICDUEQBed3', ...
        };

gfiles = {'/temp/ulbl_p/AUG/SHOTS/33353/g33353.2900_EQH',...
          '/temp/ulbl_p/AUG/WLS/33353/g33353.2900_AUGD_IDE_ed0',...
          '/temp/ulbl_p/AUG/WLS/33353/g33353.2900_MICDU_EQB_ed0',...
          '/temp/ulbl_p/AUG/WLS/33353/g33353.2900_MICDU_EQB_ed1',...
          '/temp/ulbl_p/AUG/WLS/33353/g33353.2900_MICDU_EQB_ed2',...
          %'/temp/ulbl_p/AUG/WLS/33353/g33353.2900_MICDU_EQB_ed3',...
          };

m = 4:9;
n = 2 .* ones(size(m));

shot = 33353;
time = 2900;

for k = 1:numel(equis)

    runpath = ['/temp/ulbl_p/BALANCE_2020/', [studyname, equis{k}], '/', num2str(shot), '_', num2str(time),'/'];
    
    %input
    gfile  = gfiles{k};
    filehead = ['/temp/ulbl_p/AUG/SHOTS/',num2str(shot),'/',num2str(shot),'.',num2str(time)];
    cfile  = [filehead,'_coil.dat'];

    dapath = '/temp/ulbl_p/DA/ASTRA/';
    
    neprof = [filehead,'_ne_PED_ULBLP_rho_pol.dat'];
    Teprof = [filehead,'_Te_PED_ULBLP_rho_pol.dat'];
    Tiprof = [filehead,'_Ti_PED_ULBLP_rho_pol.dat'];
    vtprof = [filehead,'_vt_PED_ULBLP_rho_pol.dat'];
    
    %location of field.dat
    pfile = ['/temp/ulbl_p/MESH3D/',num2str(shot),'/field.dat']; %will be calculated if not present
    fluxdatapath = ['/temp/ulbl_p/FLUXDATA/',equis{k},'/',num2str(time),'/']; %will be calculated if not present

    %##########################################################################
    % RUN BALANCE WITH MATLAB CLASS
    %##########################################################################
    
    gpecpath = ['/temp/ulbl_p/GPEC/', studyname, '/', equis{k}, '/', num2str(shot), '_', num2str(time),'/'];
    gpecfile = [gpecpath, 'gpec_profile_output_n2.nc'];
    
    bal = Balance(runpath, shot, time, studyname);
    bal.FLAG_FORCE_TMHD = FLAG_FORCE_TMHD;
    bal.setModes(m, n);
    bal.setCoil(pfile, cfile);
    bal.setEqui(gfile, fluxdatapath);
    bal.setTMHDCode('GPEC', gpecpath);
    bal.setProfiles(neprof, Teprof, Tiprof, vtprof);
    bal.setKiLCA();
    bal.setDaEstimation(dapath);
    bal.setOptions();
    bal.write();
    bal.run();

    %##########################################################################
    % EXPORT 2 HDF5 FORMAT
    %##########################################################################

    system(['mkdir -p ~/Balance_Results/', studyname, '/']);
    bal.export2HDF5(['~/Balance_Results/', studyname, '/'], [studyname, equis{k}]);
end

cd(mpath)