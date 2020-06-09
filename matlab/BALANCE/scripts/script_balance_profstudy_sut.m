%##########################################################################
% script_balance_equistudy.m
%##########################################################################
% description:
%--------------------------------------------------------------------------
% this script runs the balance code for the new profile study where we want
% to see the influence of the form of the n,Te,Ti,vt profiles.
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

studyname = 'ProfStudySut';
system(['mkdir -p ~/Balance_Results/', studyname, '/']);

runs = {'33353_2900_ULBLP', '33353_2900_SUT',...
        '34381_2400_ULBLP', '34381_2400_SUT',...
        '34835_4975_ULBLP', '34835_4975_SUT',...
       };

m = 4:9;
n = 2 .* ones(size(m));

for k = 1:numel(runs)

    run = strsplit(runs{k}, '_');
    shot = str2double(run{1});
    time = str2double(run{2});

    runpath = ['/temp/ulbl_p/BALANCE_2020/', studyname, '/',runs{k},'/', num2str(shot), '_', num2str(time),'/'];
    
    %input
    gfile  = ['/temp/ulbl_p/AUG/SHOTS/',num2str(shot),'/g',num2str(shot),'.',num2str(time),'_EQH'];
    cfile  = ['/temp/ulbl_p/AUG/COIL/', num2str(shot),'/',num2str(shot),'.',num2str(time),'_coil_PSL.dat'];
    dapath = '/temp/ulbl_p/DA/ASTRA/';
    
    if(strcmp(run{3}, 'ULBLP'))
        filehead = ['/temp/ulbl_p/AUG/SHOTS/',num2str(shot),'/',num2str(shot),'.',num2str(time)];
        neprof = [filehead,'_ne_PED_ULBLP_rho_pol.dat'];
        Teprof = [filehead,'_Te_PED_ULBLP_rho_pol.dat'];
        Tiprof = [filehead,'_Ti_PED_ULBLP_rho_pol.dat'];
        vtprof = [filehead,'_vt_PED_ULBLP_rho_pol.dat'];
    else
        filehead = ['/temp/ulbl_p/AUG/Suttrop/Shots/', num2str(shot), '/', num2str(time), '/',num2str(shot),'.',num2str(time)];
        neprof = [filehead,'_ne_PED.dat'];
        Teprof = [filehead,'_Te_PED.dat'];
        Tiprof = [filehead,'_Ti_PED.dat'];
        vtprof = [filehead,'_vt_PED.dat'];
    end
    %location of field.dat
    pfile = ['/temp/ulbl_p/MESH3D/',num2str(shot),'/field.dat']; %will be calculated if not present
    fluxdatapath = ['/temp/ulbl_p/FLUXDATA/',num2str(shot),'/',num2str(time),'/']; %will be calculated if not present

    %##########################################################################
    % RUN BALANCE WITH MATLAB CLASS
    %##########################################################################

    %BALANCE CODE
    bal = Balance(runpath, shot, time, [studyname, '_', runs{k}]);
    bal.FLAG_FORCE_PROFILES = FLAG_FORCE_PROFILES;
    bal.FLAG_FORCE_NEO2 = FLAG_FORCE_NEO2;
    bal.setCoil(cfile, pfile);
    bal.setEqui(gfile, fluxdatapath);
    bal.setProfiles(neprof, Teprof, Tiprof, vtprof);
    bal.setDaEstimation(dapath);
    bal.setModes(m, n);
    bal.write();
    bal.run();

    %##########################################################################
    % RUN GPEC
    %##########################################################################

    gpecpath = ['/temp/ulbl_p/GPEC/BifurcThres/', num2str(shot), '_', num2str(time), '/'];
    gpecfile = [gpecpath, 'gpec_profile_output_n2.nc'];
    
    if(~exist(gpecfile,'file') || FLAG_FORCE_GPEC)

        gpec = GPEC_interface(gpecpath, shot, time, [studyname, '_', runs{k}]);
        gpec.setEqui(gfile);
        gpec.setCoil(cfile);
        gpec.write();

        try
            gpec.run();
        catch
            warning('error in GPEC occurred.')
            gpec.run();
        end
    end
    
    %##########################################################################
    % RESCALE
    %##########################################################################

    bal.rescaleI();
    
    bal.rescaleD(gpecfile, 'GPEC');

    %##########################################################################
    % EXPORT 2 HDF5 FORMAT
    %##########################################################################

    system(['mkdir -p ~/Balance_Results/', studyname, '/']);
    bal.export2HDF5(['~/Balance_Results/', studyname, '/'], [studyname, '_', runs{k}]);
end
