%##########################################################################
% script_balance.m
%##########################################################################
% description:
%--------------------------------------------------------------------------
% this script runs the balance code for the n and Te parameterstudy with
% suttrops profiles.
%##########################################################################

%author:   Philipp Ulbl
%created:  08.05.2020

libKiLCA = '~/KiLCA_interface/';
libBalance = '~/BALANCE/balance';

addpath(genpath(libKiLCA))
addpath(genpath(libBalance))

mpath = pwd();

studyname = 'nTStudySut';
system(['mkdir -p ~/Balance_Results/', studyname, '/']);

%Runs to make
shot = 33353;
time = 2900;

m = 4:9;
n = 2 .* ones(size(m));

%##########################################################################
% 1) STANDARD RUN
%##########################################################################

runname = 'Tfac1_nfac1';
runpath = ['/temp/ulbl_p/BALANCE_2020/', studyname,'/', num2str(shot), '_', num2str(time),'/',runname,'/'];

%input
gfile  = ['/temp/ulbl_p/AUG/SHOTS/',num2str(shot),'/g',num2str(shot),'.',num2str(time),'_EQH'];
cfile  = ['/temp/ulbl_p/AUG/SHOTS/',num2str(shot),'/',num2str(shot),'.',num2str(time),'_coil.dat'];
dapath = '/temp/ulbl_p/DA/ASTRA/';

filehead = ['/temp/ulbl_p/AUG/Suttrop/Shots/',num2str(shot),'/',num2str(time),'/',num2str(shot),'.',num2str(time)];
neprof = [filehead,'_ne_PED.dat'];
Teprof = [filehead,'_Te_PED.dat'];
Tiprof = [filehead,'_Ti_PED.dat'];
vtprof = [filehead,'_vt_PED.dat'];

%location of field.dat
pfile = ['/temp/ulbl_p/MESH3D/',num2str(shot),'/field.dat']; %will be calculated if not present
fluxdatapath = ['/temp/ulbl_p/FLUXDATA/',num2str(shot),'/',num2str(time),'/']; %will be calculated if not present

%BALANCE CODE
bal = Balance(runpath, shot, time, [studyname, '_', runname]);
bal.FLAG_FORCE_PROFILES = true;
bal.setCoil(cfile, pfile);
bal.setEqui(gfile, fluxdatapath);
bal.setProfiles(neprof, Teprof, Tiprof, vtprof);
bal.setDaEstimation(dapath);
bal.setModes(m, n);
bal.write();

bal.run();

%rescale
gpecpath = ['/temp/ulbl_p/GPEC/BifurcThres/', num2str(shot), '_', num2str(time), '/'];
gpecfile = [gpecpath, 'gpec_profile_output_n2.nc'];
bal.rescaleI();
bal.rescaleD(gpecfile, 'GPEC');

system(['mkdir -p ~/Balance_Results/', studyname, '/']);
bal.export2HDF5(['~/Balance_Results/', studyname, '/ref/'], [studyname, '_', runname]);
        
%##########################################################################
% 2) VARY PROFILES
%##########################################################################

refpath = runpath;

% nfac = linspace(0.1,3,21);
% Tfac = linspace(0.5,2,31);
nfac = linspace(0,3,10);
Tfac = linspace(0,2,9);
nfac(nfac==0)=[];
Tfac(Tfac==0)=[];

for o = 1:numel(nfac)

    for p = 1:numel(Tfac)
        
        runname = ['Tfac', num2str(Tfac(p)), '_nfac', num2str(nfac(o))];
        runpath = ['/temp/ulbl_p/BALANCE_2020/', studyname,'/', num2str(shot), '_', num2str(time),'/',runname,'/'];
    
        system(['mkdir -p ', runpath]);
        system(['rsync -a ', refpath, 'profiles/ ', runpath, 'profiles/']);

        numrun = ['(',num2str((o-1)*numel(Tfac)+p),'/',num2str(numel(nfac)*numel(Tfac)),') '];
        disp([numrun, studyname, ' nfac = ', num2str(nfac(o)), ' Tfac = ', num2str(Tfac(p))]);

        %BALANCE CODE
        bal = Balance(runpath, shot, time, [studyname, '_', runname]);
        bal.setCoil(pfile, cfile);
        bal.setEqui(gfile, fluxdatapath);
        bal.setProfiles(neprof, Teprof, Tiprof, vtprof);
        bal.setModes(m, n);
        bal.write();
        
        %change density profile
        ne = load([refpath, 'profiles/n.dat']);
        ne(:, 2) = ne(:, 2) .* nfac(o);
        fid = fopen([runpath, 'profiles/n.dat'],'w');
        fprintf(fid, '%.15e\t%.15e\n', ne');
        fclose(fid);

        %change temp profile
        Te = load([refpath, 'profiles/Te.dat']);
        Te(:, 2) = Te(:, 2) .* Tfac(p);
        fid = fopen([runpath, 'profiles/Te.dat'],'w');
        fprintf(fid, '%.15e\t%.15e\n', Te');
        fclose(fid);
        
        bal.setDaEstimation(dapath);
        
        bal.run();

        %##########################################################################
        % RESCALE
        %##########################################################################

        gpecpath = ['/temp/ulbl_p/GPEC/BifurcThres/', num2str(shot), '_', num2str(time), '/'];
        gpecfile = [gpecpath, 'gpec_profile_output_n2.nc'];

        bal.rescaleI();

        bal.rescaleD(gpecfile, 'GPEC');

        %##########################################################################
        % EXPORT 2 HDF5 FORMAT
        %##########################################################################

        system(['mkdir -p ~/Balance_Results/', studyname, '/']);
        bal.export2HDF5(['~/Balance_Results/', studyname, '/'], [studyname, '_', runname]);
    end
end

cd(mpath)