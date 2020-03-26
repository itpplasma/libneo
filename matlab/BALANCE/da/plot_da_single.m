%##########################################################################
% estimate_da.m
%##########################################################################
% description:
%--------------------------------------------------------------------------
% 
%##########################################################################

addpath('~/GITHUB/libneo/matlab/Utility/xxaxis/')

%author:   Philipp Ulbl
%created:  27.02.2020

shot = 33353;
time = 2900;

%load equil file for psi pol axis
fluxdatapath = '/temp/ulbl_p/FLUXDATA/33353/2900/';
equilrqpsi = dlmread([fluxdatapath, 'equil_r_q_psi.dat'], '', 3, 0);
r = equilrqpsi(:, 1); %equivalent radius
q = equilrqpsi(:, 2); %safety factor
psi_pol_norm = equilrqpsi(:, 3)./equilrqpsi(end, 3); %normalized poloidal flux

figure('units', 'normalized', 'outerposition', [0, 0, 1, 0.7]);
axis tight

dat = load(['./out/', num2str(shot), '_', num2str(time), '_Da.dat']);

semilogy(dat(:, 1), dat(:, 2), '-', 'LineWidth', 2, ...
    'DisplayName', [num2str(shot), ' @', num2str(time), 'ms']);
hold on

xlim([40, 65])
semilogy(xlim, 1e4 .* [1, 1], '--k', 'LineWidth', 2, 'DisplayName', '"universal constant" 10^4')
ylim([1e3, 1e6])
legend()
xlabel('r / cm')
ylabel('D_a / cm^2 s^{-1}')
title('Anomalous Diffusion Coefficient Estimation')

rrho2axis(gca, xxaxis, r, psi_pol_norm, 0.06);

set(findall(gcf,'-property','FontSize'), 'FontSize', 18)
    
print('./Da_estimation_single.png', '-dpng', '-r200')
print('./Da_estimation_single.svg', '-dsvg')