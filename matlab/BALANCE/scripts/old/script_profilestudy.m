%##########################################################################
% script_profilestudy.m
%##########################################################################
% description:
%--------------------------------------------------------------------------
% this script runs the balancecode for different configurations of input
% profiles. legacy.
%##########################################################################
     
%author:   Philipp Ulbl
%created:  24.01.2020

libKiLCA = '~/KiLCA_interface/';
libBalance = '~/BALANCE/';

addpath(libKiLCA)
addpath(libBalance)
addpath([libBalance, 'balance/'])

mpath = pwd();

studyname = 'ProfileStudy_33353';
runname = 'New_Ti_vt';

 neprof = '/proj/plasma/RMP/DATA2017/33353/2.9s/neprof_aug_33353_2.9.asc';
%neprof = '/temp/ulbl_p/AUG/SHOTS/33353/33353.2900_ne_PED_ulbl_rho_pol.dat';
 Teprof = '/proj/plasma/RMP/DATA2017/33353/2.9s/Teprof_aug_33353_2.9.asc';
%Teprof = '/temp/ulbl_p/AUG/SHOTS/33353/33353.2900_Te_PED_ulbl_rho_pol.dat';
% Tiprof = '/proj/plasma/RMP/DATA2017/33353/2.9s/Tiprof_aug_33353_2.9.asc';
Tiprof = '/temp/ulbl_p/AUG/SHOTS/33353/33353.2900_Ti_PED_ulbl_rho_pol.dat';
% vtprof = '/proj/plasma/RMP/DATA2017/33353/2.9s/vtprof_aug_33353_2.9.asc';
vtprof = '/temp/ulbl_p/AUG/SHOTS/33353/33353.2900_vt_PED_ulbl_rho_pol.dat';

system(['mkdir -p ~/Balance_Results/', studyname, '/']);

%##########################################################################
% SHOT PARAMETERS
%##########################################################################

shot = 33353;
time = 2900;

runpath = ['/temp/ulbl_p/BALANCE_2020/', studyname, '/', runname, '/', num2str(shot), '_', num2str(time),'/'];

m = 4:9;
n = 2 .* ones(size(m));

%experimental input
gfile  = ['/temp/ulbl_p/AUG/SHOTS/',num2str(shot),'/g',num2str(shot),'.',num2str(time),'_EQH'];
pfile = ['/temp/ulbl_p/MESH3D/',num2str(shot),'/field.dat']; %will be calculated if not present
fluxdatapath = ['/temp/ulbl_p/FLUXDATA/',num2str(shot),'/',num2str(time),'/']; %will be calculated if not present

%##########################################################################
% RUN BALANCE WITH MATLAB CLASS
%##########################################################################

% system(['mkdir -p ', runpath]);
% system(['rsync -av ', runprof, ' ', runpath, 'profiles/']);

%BALANCE CODE
bal = Balance(runpath, shot, time, libKiLCA);
bal.FLAG_FORCE_PROFILES = true;
bal.setCoil(pfile);
bal.setEqui(gfile, fluxdatapath);
bal.setProfiles(neprof, Teprof, Tiprof, vtprof);
bal.setModes(m, n);
bal.write();
bal.run();

%##########################################################################
% BALANCE CURRENTS
%##########################################################################

Brvac = zeros(size(m));
for l=1:numel(m)
    post = bal.kil_vac.postprocessors{l};
    Brvac(l) = interp1(post.r, post.Br, post.rres);
end
A_ref = (bal.r_res .* bal.r_big ./ bal.n) .* Brvac;

bal.readAntenna();
A_realistic = zeros(size(m));
for l=1:numel(m)
    A_realistic(l) = interp1(bal.Atheta_antenna{l}(:, 1), bal.Atheta_antenna{l}(:, 2), bal.r_res(l));
end

bal.rescaleI(abs(A_realistic ./ A_ref));

%##########################################################################
% READ GPEC CURRENTS
%##########################################################################

gpecpath = ['/temp/ulbl_p/GPEC/', num2str(shot), '_', num2str(time), '/'];

qraw = ncread([gpecpath, 'gpec_profile_output_n2.nc'], 'q_rational'); %read I res in A
Iraw = ncread([gpecpath, 'gpec_profile_output_n2.nc'], 'I_res'); %read I res in A
Igpec = Iraw(:, 1) + 1i .* Iraw(:, 2);
ind = ismember(abs(qraw), abs(m./n));
Igpec = abs(Igpec(ind)') / 10; %cgs

bal.rescaleD((Igpec./bal.I_res_resc).^2);

%##########################################################################
% EXPORT 2 HDF5 FORMAT
%##########################################################################

bal.export2HDF5(['~/Balance_Results/', studyname, '/'], [studyname, '_', runname]);

%##########################################################################
% RUN REFERENCE RUN
%##########################################################################
% 
% gfile = '/proj/plasma/RMP/DATA2017/33133/3.0s/g33133.3000_ed4';
% 
% refrunpath = [runpath, 'reference/'];
% system(['mkdir -p ', refrunpath]);
% system(['rsync -av ', refprof, ' ', refrunpath, 'profiles/']);
% 
% bal_ref = Balance(refrunpath, shot, time, libKiLCA);
% bal_ref.FLAG_FORCE_PROFILES = false;
% bal_ref.setCoil(pfile);
% bal_ref.setEqui(gfile, fluxdatapath);
% bal_ref.setProfiles();
% bal_ref.setModes(m, n);
% bal_ref.write();
% bal_ref.run();
% 
% %##########################################################################
% % PLOT
% %##########################################################################
% 
% plotrefname = strrep(refname, '_', ' ');
% plotname = strrep(runname, '_', ' ');
% 
% figure('units', 'normalized', 'outerposition', [0, 0, 1, 1]);
% axis tight
% bal.plotDe22Single('DisplayName', plotname);
% bal_ref.plotDe22Single('DisplayName', plotrefname);
% legend('Location', 'northoutside')
% 
% print(['~/Balance_Results/profile_study/', filename, '_De22'], '-dpng', '-r200')
% print(['~/Balance_Results/profile_study/', filename, '_De22'], '-dsvg')
% 
% figure('units', 'normalized', 'outerposition', [0, 0, 1, 1]);
% axis tight
% bal.plotProfiles('Color', 'r', 'Handlevisibility', 'on');
% bal_ref.plotProfiles('Color', 'b', 'Handlevisibility', 'on');
% ah = gca;
% for k=(numel(ah.Children)-1):-1:(numel(ah.Children)/2+1)
%     ah.Children(k).HandleVisibility = 'off';
% end
% s = arrayfun(@(i) ['m=', num2str(m(i)), ', n=', num2str(n(i))], 1:numel(m), 'UniformOutput', false);
% legend(plotname, plotrefname, s{:})
% 
% print(['~/Balance_Results/profile_study/', filename, '_Prof'], '-dpng', '-r200')
% print(['~/Balance_Results/profile_study/', filename, '_Prof'], '-dsvg')







