%##########################################################################
% script_balance.m
%##########################################################################
% description:
%--------------------------------------------------------------------------
% this script runs the balance code for the n and Te parameterstudy with
% suttrops profiles.
%##########################################################################

%author:   Philipp Ulbl
%created:  15.06.2020

libKiLCA = '~/KiLCA_interface/';
libBalance = '~/BALANCE/balance';

addpath(genpath(libKiLCA))
addpath(genpath(libBalance))

mpath = pwd();

studyname = 'VelScan';
system(['mkdir -p ~/Balance_Results/', studyname, '/']);

%Runs to make
shot = 33133;
time = 3000;

m = 5:7;
n = 2 .* ones(size(m));

%##########################################################################
% 1) STANDARD RUN
%##########################################################################

runname = 'Vzfac1';
runpath = ['/temp/ulbl_p/BALANCE_2020/', studyname,'/', num2str(shot), '_', num2str(time),'/',runname,'/'];
refpath = runpath;

%input
gfile  = ['/temp/ulbl_p/AUG/SHOTS/',num2str(shot),'/g',num2str(shot),'.',num2str(time),'_EQH'];
cfile  = ['/temp/ulbl_p/AUG/SHOTS/',num2str(shot),'/',num2str(shot),'.',num2str(time),'_coil.dat'];
dapath = '/temp/ulbl_p/DA/ASTRA/';

filehead = ['/temp/ulbl_p/AUG/SHOTS/',num2str(shot),'/',num2str(shot),'.',num2str(time)];
neprof = [filehead,'_ne_PED_ULBLP_rho_pol.dat'];
Teprof = [filehead,'_Te_PED_ULBLP_rho_pol.dat'];
Tiprof = [filehead,'_Ti_PED_ULBLP_rho_pol.dat'];
vtprof = [filehead,'_vt_PED_ULBLP_rho_pol.dat'];

%location of field.dat
fluxdatapath = ['/temp/ulbl_p/FLUXDATA/',num2str(shot),'/',num2str(time),'/']; %will be calculated if not present

gpecpath = ['/temp/ulbl_p/GPEC/TimeEvol/', num2str(shot), '_', num2str(time),'/'];

%BALANCE CODE
bal = Balance(runpath, shot, time, [studyname, '_', runname]);
bal.setModes(m, n);
bal.setCoil(cfile);
bal.setEqui(gfile, fluxdatapath);
bal.setTMHDCode('GPEC', gpecpath);
bal.setProfiles(neprof, Teprof, Tiprof, vtprof);
bal.setKiLCA();
bal.setDaEstimation(dapath);
bal.write();
bal.run();

system(['mkdir -p ~/Balance_Results/', studyname, '/ref/']);
bal.export2HDF5(['~/Balance_Results/', studyname, '/ref/'], [studyname, '_', runname]);
        
%##########################################################################
% 2) VARY PROFILES
%##########################################################################

Vzfac = linspace(-0.5, 10, 526);
Vzfac=unique(Vzfac);

for o = 366:numel(Vzfac)

    runname = ['Vzfac', sprintf('%.2f',Vzfac(o))];
    runpath = ['/temp/ulbl_p/BALANCE_2020/', studyname,'/', num2str(shot), '_', num2str(time),'/',runname,'/'];

    system(['mkdir -p ', runpath]);
    system(['rsync -a ', refpath, 'profiles/ ', runpath, 'profiles/']);

    numrun = ['(',num2str(o),'/',num2str(numel(Vzfac)),') '];
    disp([numrun, studyname, ' Vzfac = ', num2str(Vzfac(o)),]);

    %BALANCE CODE
    bal = Balance(runpath, shot, time, [studyname, '_', runname]);
    bal.setModes(m, n);
    bal.setCoil(cfile);
    bal.setEqui(gfile, fluxdatapath);
    bal.setTMHDCode('GPEC', gpecpath);
    bal.setProfiles(neprof, Teprof, Tiprof, vtprof);
    bal.setKiLCA();
    bal.write();

    %change temp profile
    Vz = load([refpath, 'profiles/Vz.dat']);
    Vz(:, 2) = Vz(:, 2) .* Vzfac(o);
    fid = fopen([runpath, 'profiles/Vz.dat'],'w');
    fprintf(fid, '%.15e\t%.15e\n', Vz');
    fclose(fid);

    bal.setDaEstimation(dapath);
    bal.run();

    %##########################################################################
    % EXPORT 2 HDF5 FORMAT
    %##########################################################################

    system(['mkdir -p ~/Balance_Results/', studyname, '/']);
    bal.export2HDF5(['~/Balance_Results/', studyname, '/'], [studyname, '_', runname]);
    
end

cd(mpath)