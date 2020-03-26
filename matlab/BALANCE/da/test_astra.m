%##########################################################################
% estimate_da.m
%##########################################################################
% description:
%--------------------------------------------------------------------------
% 
%##########################################################################

addpath('~/GITHUB/libneo/matlab/Utility/xxaxis/')
addpath('~/GITHUB/libneo/matlab/Utility/smooth2level/')

%author:   Philipp Ulbl
%created:  19.03.2020

shot = 33353;
time = 2900;

profpath = ['/temp/ulbl_p/BALANCE_2020/suppression/', ...
    num2str(shot), '_', num2str(time), '/profiles/'];
    
%load equil file for psi pol axis
fluxdatapath = '/temp/ulbl_p/FLUXDATA/33353/2900/';
equilrqpsi = dlmread([fluxdatapath, 'equil_r_q_psi.dat'], '', 3, 0);
r = equilrqpsi(:, 1); %equivalent radius
q = equilrqpsi(:, 2); %safety factor
Vequi = equilrqpsi(:, 7);
psi_pol_norm = equilrqpsi(:, 3)./equilrqpsi(end, 3); %normalized poloidal flux

e = 1.602e-19;

%##########################################################################
% LOAD ASTRA DATA
%##########################################################################

raw = importdata(['/temp/ulbl_p/DA/ASTRA/ASTRA_PB_AUG_#', num2str(shot),'_t', num2str(time/1000)]);
dat = raw.data';
dat = dat(~isnan(dat));
dat = reshape(dat, 19, numel(dat)/19)';

rho_pol = dat(:, 3);
chi_e = interp1(rho_pol.^2, dat(:, 4), psi_pol_norm, 'pchip') .* 1e4;
Q_e = interp1(rho_pol.^2, dat(:, 6), psi_pol_norm, 'pchip') .* 1e6 ./ e;
T_e = interp1(rho_pol.^2, dat(:, 8), psi_pol_norm, 'pchip') .* 1e3;
n_e = interp1(rho_pol.^2, dat(:, 10), psi_pol_norm, 'pchip') .* 1e13;
V = interp1(rho_pol.^2, dat(:, 13), psi_pol_norm, 'pchip') .* 1e6;

%##########################################################################
% RECONSTRUCT DA
%##########################################################################

smo = smooth2level(T_e, r, 2, 5e-2);
T_e = smo{1};
dT_e = smo{2};

smo = smooth2level(V, r, 2, 5e-2);
V = smo{1};
S = smo{2};
Da_e = -Q_e ./ n_e ./ dT_e ./ S;

%##########################################################################
% ESTIMATE WITH OWN PROFILES
%##########################################################################

Sequi = gradient(Vequi, r);

ne = load([profpath, 'n.dat']);
Te = load([profpath, 'Te.dat']);
ne_prof = interp1(ne(:, 1), ne(:, 2), r, 'pchip');
Te_prof = interp1(Te(:, 1), Te(:, 2), r, 'pchip');
dT_prof = gradient(Te_prof, r);

Da_prof = -Q_e ./ ne_prof ./ dT_prof ./ Sequi;

%##########################################################################
% PLOT
%##########################################################################

figure('units', 'normalized', 'outerposition', [0, 0, 1, 0.7]);
axis tight

dat = load(['./out/', num2str(shot), '_', num2str(time), '_Da.dat']);
semilogy(dat(:, 1), dat(:, 2), '-', 'LineWidth', 2, ...
    'DisplayName', 'Estimation');
hold on
semilogy(r, chi_e, '-', 'LineWidth', 2, ...
    'DisplayName', 'ASTRA \chi_e');
semilogy(r, Da_e, '-', 'LineWidth', 2, ...
    'DisplayName', 'ASTRA D_a reconstructed');
semilogy(r, Da_prof, '-', 'LineWidth', 2, ...
    'DisplayName', 'Estimation with ASTRA P_{inp}=QETOT');

%xlim([40, 65])
semilogy(xlim, 1e4 .* [1, 1], '--k', 'LineWidth', 2, 'DisplayName', '"universal constant" 10^4')
ylim([1e3, 1e6])
xlabel('r / cm')
ylabel('D_a / cm^2 s^{-1}')
title('Anomalous Diffusion Coefficient')

set(gca, 'outerposition', [0, 0, 0.8, 1])
legend('Location', [0.85, 0.1, 0, 1])

rrho2axis(gca, xxaxis, r, psi_pol_norm, 0.2);

set(findall(gcf,'-property','FontSize'), 'FontSize', 18)
    
print('./Da_compare_ASTRA.png', '-dpng', '-r200')