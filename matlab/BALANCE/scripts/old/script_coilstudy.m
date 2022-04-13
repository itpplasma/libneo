%##########################################################################
% script_coilstudy.m
%##########################################################################
% description:
%--------------------------------------------------------------------------
% this script runs the balance code for a legacy run with old and new coil
% file and prints the difference.
%##########################################################################

%author:   Philipp Ulbl
%created:  17.01.2020

libKiLCA = '~/KiLCA_interface/';
libBalance = '~/BALANCE/';

addpath(libKiLCA)
addpath(libBalance)
addpath([libBalance, 'balance/'])

mpath = pwd();

m = 4:10;
n = 2 .* ones(size(m));

%##########################################################################
% AUG PARAMETERS
%##########################################################################

%these values are a bit arbitrary but must be larger than r_end
%r_end is the maximum radius in equil_r_q_psi.dat which will be calculated
%below ( = real position of separatrix)
r_sep = 67;             %position of separatrix
r_ant = 70;             %position of antenna
r_idw = 80;             %position of ideal wall

r_sta = 3;              %minimum radius in calculations of KiLCA
r_min = 0.1;            %minimum radius to take from exp data

r_max = r_ant;          %maximum radius in balance
r_max_prof = r_max+1;   %maximum radius for profiles (> r_ant)

rb_cut_out = 67.5;      %needed in balance.in
re_cut_out = 68;        %needed in balance.in

%location of convexfile for AUG machine
convexfile = '/proj/plasma/RMP/DATA/convexwall.dat';
%location of math libraries
softpath = '/proj/plasma/soft/';

%##########################################################################
% LOAD EXISTING DATA
%##########################################################################

shot = 33133;
time = 2600;

runpath = ['/temp/ulbl_p/BALANCE_2017/CoilStudy_Remake33133_', num2str(shot), '_', num2str(time),'/'];
profpath = [runpath, 'profiles/'];

%if these 2 lines not work you have to specify location by hand
time_s = sprintf('%1.1f', time/1000);
fluxdatapath_sergei = ['/proj/plasma/RMP/DATA2017/', num2str(shot),'/',time_s,'s/PREPROCESS/'];
fluxdatapath_philipp = '/temp/ulbl_p/FLUXDATA/RemakeLegacy/33133/2600/';
profload = ['/temp/heyn/BALANCE_2017/',num2str(shot),'_',num2str(time),'/profiles/'];

[~, ~] = system(['rsync -av ', profload, ' ', profpath]);
[~, ~] = system(['rm ', profpath, 'Vz.dat']);
[~, ~] = system(['cp ', profpath, 'Vz_org.dat ', profpath, 'Vz.dat']);

%##########################################################################
% CREATE RUN PATH
%##########################################################################

system(['mkdir -p ', runpath]);
[~, ~] = system(['rsync -av ', libBalance, 'template/ ', runpath]);
%##########################################################################
% FOURIERMODES
%##########################################################################

%needed output files
amnpath = [fluxdatapath_sergei, 'amn.dat'];
btorrbigpath = [fluxdatapath_sergei, 'btor_rbig.dat'];
equilrqpsipath = [fluxdatapath_sergei, 'equil_r_q_psi.dat'];

%output of fouriermodes: rbig, btor
raw = load(btorrbigpath);
b_tor = raw(1);
r_big = raw(2);

%output of fouriermodes: equil r q psi
%read equil file
raw = equil_read(equilrqpsipath);
%extract r, q, psi
r = raw(:,1);                           % equivalent radius
q = raw(:,2);                           % safety factor q
psi_pol_norm = raw(:,3)./raw(end,3);    % normalized poloidal flux

%##########################################################################
% PREPARE KiLCA - Sergei
%##########################################################################

%FLRE, VACUUM USING KiLCA CLASS
nmodes = 1;
rad = [r_sta, r_sep, r_ant, r_idw];
bound = {'center', 'interface', 'antenna', 'idealwall'};
med = {'vacuum', 'vacuum', 'vacuum'};

kil_vac_sergei = KiLCA_interface([runpath, 'kilca_sergei/'], 'vacuum');
kil_vac_sergei.BLUE_PATH = [libKiLCA, 'blueprints/'];
kil_vac_sergei.PROF_PATH = profpath;
kil_vac_sergei.set_background(r_big, r_sep);
kil_vac_sergei.background.Btor = b_tor;
kil_vac_sergei.set_antenna(r_ant, nmodes);
med{1} = kil_vac_sergei.run_type;
kil_vac_sergei.set_zones(rad, bound, med);

kil_flre_sergei = KiLCA_interface([runpath, 'kilca_sergei/'], 'flre');
kil_flre_sergei.BLUE_PATH = [libKiLCA, 'blueprints/'];
kil_flre_sergei.PROF_PATH = profpath;
kil_flre_sergei.set_background(r_big, r_sep);
kil_flre_sergei.background.Btor = b_tor;
kil_flre_sergei.set_antenna(r_ant, nmodes);
med{1} = kil_flre_sergei.run_type;
kil_flre_sergei.set_zones(rad, bound, med);

kil_vac_sergei.antenna.I0 = 1e13; %like in Martins runs
kil_vac_sergei.antenna.flab(1) = 1; %1Hz
kil_vac_sergei.zones{3}.vacuum.sigma(1) = 1.3e16; %resistive wall

kil_flre_sergei.antenna.I0 = 1e13; %like in Martins runs
kil_flre_sergei.antenna.flab(1) = 1; %1Hz
kil_flre_sergei.zones{3}.vacuum.sigma(1) = 1.3e16; %resistive wall

%##########################################################################
% BALANCE OPTIONS
%##########################################################################

%BALANCE OPTIONS
balopt_sergei = balanceoptions(kil_flre_sergei.pathofrun, kil_vac_sergei.pathofrun);
balopt_sergei.Btor = b_tor;
balopt_sergei.Rtor = r_big;
balopt_sergei.rmin = r_sta;
balopt_sergei.rmax = r_ant;
balopt_sergei.rsep = r_sep;
balopt_sergei.rb_cut_out = rb_cut_out;
balopt_sergei.re_cut_out = re_cut_out;

%##########################################################################
% RUN BALANCE WITH FLUXDATA BY SERGEI
%##########################################################################

system(['ln -sfT ', kil_vac_sergei.path, 'vacuum/ ', runpath, 'vacuum']);
system(['ln -sfT ', kil_flre_sergei.path, 'flre/ ', runpath, 'flre']);

%FIELD DIVB0
fdb0 = field_divB0('', '', '', fluxdatapath_sergei);
fdb0.write([libBalance, 'blueprints/'], runpath);
system(['ln -sfT ', softpath, ' ', runpath, 'soft']);

%BALANCE CODE
bal_sergei = balancecode(runpath, softpath);
bal_sergei.BLUE_PATH = [libBalance, 'blueprints/'];
bal_sergei.setOptions(balopt_sergei);
bal_sergei.setFactors(0);
bal_sergei.setModes(m, n);
bal_sergei.setKiLCA(kil_flre_sergei, kil_vac_sergei);
bal_sergei.path_output = [runpath, 'out_sergei/'];
bal_sergei.write();
bal_sergei.run();
bal_sergei.loadOutput();

%##########################################################################
% FOURIERMODES
%##########################################################################

%needed output files
amnpath = [fluxdatapath_philipp, 'amn.dat'];
btorrbigpath = [fluxdatapath_philipp, 'btor_rbig.dat'];
equilrqpsipath = [fluxdatapath_philipp, 'equil_r_q_psi.dat'];

%output of fouriermodes: rbig, btor
raw = load(btorrbigpath);
b_tor = raw(1);
r_big = raw(2);

%output of fouriermodes: equil r q psi
%read equil file
raw = equil_read(equilrqpsipath);
%extract r, q, psi
r = raw(:,1);                           % equivalent radius
q = raw(:,2);                           % safety factor q
psi_pol_norm = raw(:,3)./raw(end,3);    % normalized poloidal flux

%##########################################################################
% PREPARE KiLCA - Philipp
%##########################################################################

kil_vac_philipp = KiLCA_interface([runpath, 'kilca_philipp/'], 'vacuum');
kil_vac_philipp.BLUE_PATH = [libKiLCA, 'blueprints/'];
kil_vac_philipp.PROF_PATH = profpath;
kil_vac_philipp.set_background(r_big, r_sep);
kil_vac_philipp.background.Btor = b_tor;
kil_vac_philipp.set_antenna(r_ant, nmodes);
med{1} = kil_vac_philipp.run_type;
kil_vac_philipp.set_zones(rad, bound, med);

kil_flre_philipp = KiLCA_interface([runpath, 'kilca_philipp/'], 'flre');
kil_flre_philipp.BLUE_PATH = [libKiLCA, 'blueprints/'];
kil_flre_philipp.PROF_PATH = profpath;
kil_flre_philipp.set_background(r_big, r_sep);
kil_flre_philipp.background.Btor = b_tor;
kil_flre_philipp.set_antenna(r_ant, nmodes);
med{1} = kil_flre_philipp.run_type;
kil_flre_philipp.set_zones(rad, bound, med);

kil_vac_philipp.antenna.I0 = 1e13; %like in Martins runs
kil_vac_philipp.antenna.flab(1) = 1; %1Hz
kil_vac_philipp.zones{3}.vacuum.sigma(1) = 1.3e16; %resistive wall

kil_flre_philipp.antenna.I0 = 1e13; %like in Martins runs
kil_flre_philipp.antenna.flab(1) = 1; %1Hz
kil_flre_philipp.zones{3}.vacuum.sigma(1) = 1.3e16; %resistive wall

%##########################################################################
% BALANCE OPTIONS
%##########################################################################

%BALANCE OPTIONS
balopt_philipp = balanceoptions(kil_flre_philipp.pathofrun, kil_vac_philipp.pathofrun);
balopt_philipp.Btor = b_tor;
balopt_philipp.Rtor = r_big;
balopt_philipp.rmin = r_sta;
balopt_philipp.rmax = r_ant;
balopt_philipp.rsep = r_sep;
balopt_philipp.rb_cut_out = rb_cut_out;
balopt_philipp.re_cut_out = re_cut_out;

%##########################################################################
% RUN BALANCE WITH FLUXDATA BY PHILIPP
%##########################################################################

system(['ln -sfT ', kil_vac_philipp.path, 'vacuum/ ', runpath, 'vacuum']);
system(['ln -sfT ', kil_flre_philipp.path, 'flre/ ', runpath, 'flre']);

%FIELD DIVB0
fdb0 = field_divB0('', '', '', fluxdatapath_philipp);
fdb0.write([libBalance, 'blueprints/'], runpath);
system(['ln -sfT ', softpath, ' ', runpath, 'soft']);

%BALANCE CODE
bal_philipp = balancecode(runpath, softpath);
bal_philipp.BLUE_PATH = [libBalance, 'blueprints/'];
bal_philipp.setOptions(balopt_philipp);
bal_philipp.setFactors(0);
bal_philipp.setModes(m, n);
bal_philipp.setKiLCA(kil_flre_philipp, kil_vac_philipp);
bal_philipp.path_output = [runpath, 'out_philipp/'];
bal_philipp.write();
bal_philipp.run();
bal_philipp.loadOutput();

%##########################################################################
% POSTPROCESS
%##########################################################################

%load output of Balance
figure('units', 'normalized', 'outerposition', [0, 0, 1, 1]);
axis tight
subplot(2, 1, 1)
hold on
for i = find(m==4):find(m==5)
    f5k = bal_sergei.outputs{i}.fort5000;
    f5k.plot('de22', 'D_e22', 'LineStyle', '-', 'LineWidth', 2, 'DisplayName', ['m=', num2str(m(i)), ', n=', num2str(n(i)), ' - Sergei']);
    f5k = bal_philipp.outputs{i}.fort5000;
    f5k.plot('de22', 'D_e22', 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', ['m=', num2str(m(i)), ', n=', num2str(n(i)), ' - Philipp']);
end
plot(xlim, 1e4 .* [1, 1], '-k', 'DisplayName', 'Da')
hold off
title(['coil study: shot ', num2str(shot), ' @', num2str(time), 'ms'])
legend('Location', 'eastoutside')
xlim([46, 54])
subplot(2, 1, 2)
hold on
for i = find(m==6):find(m==9)
    f5k = bal_sergei.outputs{i}.fort5000;
    f5k.plot('de22', 'D_e22', 'LineStyle', '-', 'LineWidth', 2, 'DisplayName', ['m=', num2str(m(i)), ', n=', num2str(n(i)), ' - Sergei']);
    f5k = bal_philipp.outputs{i}.fort5000;
    f5k.plot('de22', 'D_e22', 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', ['m=', num2str(m(i)), ', n=', num2str(n(i)), ' - Philipp']);
end
plot(xlim, 1e4 .* [1, 1], '-k', 'DisplayName', 'Da')
hold off
title(['coil study: shot ', num2str(shot), ' @', num2str(time), 'ms'])
legend()
legend('Location', 'eastoutside')
xlim([56, 64])

set(findall(gcf,'-property','FontSize'), 'FontSize', 18)

%print(['~/BALANCE/coil_study/', num2str(shot), '_', num2str(time)], '-dsvg')
%print(['~/BALANCE/coil_study/', num2str(shot), '_', num2str(time)], '-dpng', '-r500')

%##########################################################################
% Electric Field
%##########################################################################

% bal_sergei.kil_flre.post();
% bal_philipp.kil_flre.post();
% 
% bal_sergei.kil_flre.multimode.plotVs(bal_philipp.kil_flre, 'Er', 'E', 'Abs');
% legend('Coils - Sergei', 'q=m/n', 'Coils - Philipp')



