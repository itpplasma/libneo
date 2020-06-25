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

clear all %IS necessary

libKiLCA = '~/KiLCA_interface/';
libBalance = '~/BALANCE/balance';

addpath(genpath(libKiLCA))
addpath(genpath(libBalance))

mpath = pwd();

studyname = 'VelScan_dTe_c5';
system(['mkdir -p ~/Balance_Results/', studyname, '/']);

%Runs to make
shot = 33133;
time = 3000;

m = 6;
n = 2 .* ones(size(m));

copy = '/temp/ulbl_p/BALANCE_2020/VelScan_dTe/33133_3000/VzfacRef/profiles/';

%##########################################################################
% 1) STANDARD RUN
%##########################################################################

runname = 'VzfacRef';
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
bal.setProfiles(neprof, Teprof, Tiprof, vtprof, copy);
bal.setKiLCA();
bal.setDaEstimation(dapath);
bal.write();
bal.run();

system(['mkdir -p ~/Balance_Results/', studyname, '/ref/']);
bal.export2HDF5(['~/Balance_Results/', studyname, '/ref/'], [studyname, '_', runname]);

%##########################################################################
% PRECOMP
%##########################################################################

%calculate factor for Er recomputation from reference run
prof = bal.profiles;

if(~isempty(prof.B_out))
    B0 = prof.B_out;
else
    B0 = interp1(bal.kil_vacuum.backgrounddata.b0(:, 1), ...
                 bal.kil_vacuum.backgrounddata.b0(:, 2), prof.r_out);
end

ErVzfac = prof.r_out .* B0 ./ (prof.r_big * prof.qp.y_out .* 3e10);

%set factors for vz scan
%Vzfac = linspace(3, 6, 101);
Vzfac = linspace(0, 6, 241);
Vzfac = unique(Vzfac);

copy = ['/temp/ulbl_p/BALANCE_2020/',studyname,'/33133_3000/VzfacRef/profiles/'];

%##########################################################################
% 2) VARY TE PROFILES
%##########################################################################

tepath = '/temp/ulbl_p/PROFILES/Te_Grad/';
dirs = dir(tepath);
dirs = dirs(3:end);
dirs = dirs([dirs.isdir]==true);

dirname = arrayfun(@(c) split(c.name,'_'), dirs, 'UniformOutput', false);
alpha = cell2mat(cellfun(@(c) str2double(c{2}), dirname, 'UniformOutput', false));

for k = 1:numel(alpha)
    
    if(alpha(k)>0.0)
        continue;
    end
    
    %##########################################################################
    % 3) VARY VZ PROFILES
    %##########################################################################

    for o = 1:numel(Vzfac)

%         if(alpha(k)==0.7 && o < 99)
%             continue;
%         end
        
        runname = ['Vzfac', sprintf('%.2f',Vzfac(o))];
        runpath = ['/temp/ulbl_p/BALANCE_2020/', studyname, '/', dirs(k).name, '/', num2str(shot), '_', num2str(time),'/',runname,'/'];

        numrun = ['(',num2str(o+(k-1)*numel(Vzfac)),'/',num2str(numel(Vzfac)*numel(alpha)),') '];
        disp([numrun, studyname, ' alpha = ', num2str(alpha(k)), ' Vzfac = ', num2str(Vzfac(o)),]);

        %BALANCE CODE
        bal.changeRun(runpath, [studyname, '_', runname]);
        bal.setProfiles(neprof, Teprof, Tiprof, vtprof, copy);
        bal.kil_vacuum.background.ce = 5;
        bal.kil_vacuum.background.ci = 5;
        bal.kil_vacuum.output.backdata = 1;
        bal.kil_vacuum.output.lindata = 1;
        bal.kil_vacuum.output.varquant = 1;
        bal.kil_flre.background.ce = 5;
        bal.kil_flre.background.ci = 5;
        bal.kil_flre.output.backdata = 1;
        bal.kil_flre.output.lindata = 1;
        bal.kil_flre.output.varquant = 1;
        bal.write();

        %change Te profile
        system(['rsync -a ', tepath, dirs(k).name, '/Te.dat ', runpath, 'profiles/Te.dat']);
        
        %change vz profile
        Vz_old = load([refpath, 'profiles/Vz.dat']);
        Vz_new = Vz_old;
        Vz_new(:, 2) = Vz_new(:, 2) .* (1 + Vzfac(o));
        fid = fopen([runpath, 'profiles/Vz.dat'],'w');
        fprintf(fid, '%.15e\t%.15e\n', Vz_new');
        fclose(fid);

        %change Er profile (recalculation with new vz)
        Er = load([refpath, 'profiles/Er.dat']);
        Er(:, 2) = Er(:, 2) + ErVzfac .* (Vz_new(:, 2) - Vz_old(:, 2));
        fid = fopen([runpath, 'profiles/Er.dat'],'w');
        fprintf(fid, '%.15e\t%.15e\n', Er');
        fclose(fid);

        %run
        bal.run();

        %##########################################################################
        % EXPORT 2 HDF5 FORMAT
        %##########################################################################

        system(['mkdir -p ~/Balance_Results/', studyname, '/', dirs(k).name, '/']);
        bal.export2HDF5(['~/Balance_Results/', studyname, '/', dirs(k).name, '/'], [studyname, '_', runname]);

    end
    
end

cd(mpath)