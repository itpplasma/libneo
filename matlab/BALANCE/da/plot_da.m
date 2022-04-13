%##########################################################################
% estimate_da.m
%##########################################################################
% description:
%--------------------------------------------------------------------------
% 
%##########################################################################
     
%author:   Philipp Ulbl
%created:  27.02.2020

%Runs to make
run = [30835, 2300;...
       30835, 3000;...
       33120, 5500;...
       33120, 5635;...
       33133, 2600;...
       33133, 3000;...
       33353, 2325;...
       33353, 2670;...
       33353, 2900];

style = {'-','-','--','--',':',':','-.','-.','-.'};
   
figure('units', 'normalized', 'outerposition', [0, 0, 1, 1]);
axis tight
for k = 1:size(run, 1)

    shot = run(k, 1);
    time = run(k, 2);
    
    dat = load(['./out/', num2str(shot), '_', num2str(time), '_Da.dat']);
    
    semilogy(dat(:, 1), dat(:, 2), style{k}, 'LineWidth', 2, ...
        'DisplayName', [num2str(run(k, 1)), ' @', num2str(run(k, 2)), 'ms']);
    hold on
end
xlim([40, max(xlim)])
semilogy(xlim, 1e4 .* [1, 1], '--k', 'LineWidth', 2, 'DisplayName', '"universal constant" 10^4')
ylim([1e1, 1e6])

legend('Location', 'southeast')
xlabel('r / cm')
ylabel('D_a / cm^2 s^{-1}')
title('Anomalous Diffusion Coefficient Estimation')

set(findall(gcf,'-property','FontSize'), 'FontSize', 18)
    
print('./Da_estimation.png', '-dpng', '-r200')
print('./Da_estimation.svg', '-dsvg')
    