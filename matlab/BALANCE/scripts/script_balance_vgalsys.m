%##########################################################################
% script_balance_vgalsys.m
%##########################################################################
% description:
%--------------------------------------------------------------------------
% this script makes a scan over different vgalsys from KiLCA for all modes.
% it outputs a plot of de22 over vgalsys for all modes
%##########################################################################

%author:   Philipp Ulbl
%created:  13.05.2020

mpath = pwd();

libKiLCA = '~/KiLCA_interface/';
libBalance = '~/BALANCE/balance';

addpath(genpath(libKiLCA))
addpath(genpath(libBalance))
addpath(genpath('~/GITHUB/libneo/matlab/GPEC_interface/'));

studyname = 'VgalsysStudy';
system(['mkdir -p ~/Balance_Results/', studyname, '/']);

m = 4:9;
n = 2 .* ones(size(m));

shot = 33120;
time = 5500;

vgalsys = -logspace(10, 5, 11);
De22 = zeros(numel(vgalsys), numel(m));

for k = 1:numel(vgalsys)

    runpath = ['/temp/ulbl_p/BALANCE_2020/', studyname, '/', num2str(shot), '_', num2str(time),'/'];
    
    %input
    gfile  = ['/temp/ulbl_p/AUG/SHOTS/',num2str(shot),'/g',num2str(shot),'.',num2str(time),'_EQH'];
    filehead = ['/temp/ulbl_p/AUG/SHOTS/',num2str(shot),'/',num2str(shot),'.',num2str(time)];
    cfile  = [filehead,'_coil.dat'];

    dapath = '/temp/ulbl_p/DA/ASTRA/';
    
    neprof = [filehead,'_ne_PED_ULBLP_rho_pol.dat'];
    Teprof = [filehead,'_Te_PED_ULBLP_rho_pol.dat'];
    Tiprof = [filehead,'_Ti_PED_ULBLP_rho_pol.dat'];
    vtprof = [filehead,'_vt_PED_ULBLP_rho_pol.dat'];
    
    %location of field.dat
    pfile = ['/temp/ulbl_p/MESH3D/',num2str(shot),'/field.dat']; %will be calculated if not present
    fluxdatapath = ['/temp/ulbl_p/FLUXDATA/',num2str(shot),'/',num2str(time),'/']; %will be calculated if not present

    %##########################################################################
    % RUN BALANCE WITH MATLAB CLASS
    %##########################################################################
    
    gpecpath = ['/temp/ulbl_p/GPEC/BifurcThres/', num2str(shot), '_', num2str(time), '/'];
    gpecfile = [gpecpath, 'gpec_profile_output_n2.nc'];
    
    bal = Balance(runpath, shot, time, studyname);
    bal.FLAG_FORCE_PROFILES = true;
    bal.KiLCA_vgalsys = vgalsys(k);
    bal.setModes(m, n);
    bal.setCoil(cfile, pfile);
    bal.setEqui(gfile, fluxdatapath);
    bal.setTMHDCode('GPEC', gpecpath);
    bal.setProfiles(neprof, Teprof, Tiprof, vtprof, '/temp/ulbl_p/BALANCE_2020/BifurcThres/33120_5500/profiles/');
    bal.setKiLCA();
    bal.setDaEstimation(dapath);
    bal.setOptions();
    bal.write();
    bal.run();
    
    De22(k, :) = bal.De22_res;
end

%##########################################################################
% EXPORT ONE OF THEM 2 HDF5 FORMAT
%##########################################################################

system(['mkdir -p ~/Balance_Results/', studyname, '/']);
bal.export2HDF5(['~/Balance_Results/', studyname, '/'], studyname);

%##########################################################################
% Make Plot
%##########################################################################
    
system(['mkdir -p ~/Balance_Results/', studyname, '/']);

figure('units', 'normalized', 'outerposition', [0, 0, 1, 1])
axis tight
for k = 1:numel(m)

    lh = loglog(vgalsys, De22(:, k), 'o-', 'DisplayName', ['m = ', num2str(m(k))]);
    set(lh, 'MarkerFaceColor', get(lh, 'Color'));
    hold on
end
legend('Location', 'SouthWest')
title(['De22(vgalsys) for shot ', num2str(shot), ' time ', num2str(time), 'ms'])
xlabel('vgalsys / cm s^{-1}')
ylabel('D_{e22}^{QL} / cm^2 s^{-1}')
set(findall(gcf,'-property','FontSize'), 'FontSize', 22)
print(['~/Balance_Results/', studyname, '/de22.png'], '-dpng', '-r200')
print(['~/Balance_Results/', studyname, '/de22.svg'], '-dsvg')
cd(mpath)
