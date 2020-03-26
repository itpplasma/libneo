%##########################################################################
% estimate_da.m
%##########################################################################
% description:
%--------------------------------------------------------------------------
% 
%##########################################################################
     
%author:   Philipp Ulbl
%created:  21.01.2020

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

%total estimated heating power from aug intranet  
PH = [30835, 2300, 8.757e6;...
      30835, 3000, 11.242e6;...
      33120, 5500, 8.675e6;...
      33120, 5635, 8.675e6;...
      33133, 2600, 9.139e6;...
      33133, 3000, 11.739e6;...
      33353, 2325, 8.568e6;...
      33353, 2670, 7.994e6;...
      33353, 2900, 9.974e6]; %in W
  
%total estimated radiated power from aug intranet
PR = [30835, 2300, 4.17e6;...
      30835, 3000, 3.54e6;...
      33120, 5500, 2.67e6;...
      33120, 5635, 2.79e6;...
      33133, 2600, 1.84e6;...
      33133, 3000, 2.61e6;...
      33353, 2325, 4.02e6;...
      33353, 2670, 2.06e6;...
      33353, 2900, 1.95e6]; %in W

for k = 1:size(run, 1)

    shot = run(k, 1);
    time = run(k, 2);

    profpath = ['/temp/ulbl_p/BALANCE_2020/suppression/', ...
        num2str(shot), '_', num2str(time), '/profiles/'];
    equilpath = ['/temp/ulbl_p/FLUXDATA/', ...
        num2str(shot), '/', num2str(time), '/equil_r_q_psi.dat'];
    btorrbigpath = ['/temp/ulbl_p/FLUXDATA/', ...
        num2str(shot), '/', num2str(time), '/btor_rbig.dat'];

    %get S from dV/dr from equil r q psi
    M = dlmread(equilpath, '', 3, 0);
    reff = M(:, 1);
    Veff = M(:, 7);
    Seff = gradient(Veff, reff);

    %##########################################################################
    % LOAD PROFILES
    %##########################################################################

    ne = load([profpath, 'n.dat']);
    Te = load([profpath, 'Te.dat']);
    r = ne(:, 1);
    ne = ne(:, 2);

    dT = gradient(Te(:, 2), Te(:, 1));

    %##########################################################################
    % CALCULATE DA
    %##########################################################################

    raw = load(btorrbigpath);            %large torus radius
    R = raw(2);
    S_est = 4*pi^2*R*r;     %flux surface area
    %spline S on r of profiles
    S = interp1(reff, Seff, r, 'linear');

    e = 1.602e-19;      %electron charge
    r_ind = PR(:,1)==shot & PR(:,2)==time;
    P_rad = PR(r_ind, 3);      %total radiated power in W
    h_ind = PH(:,1)==shot & PH(:,2)==time;
    P_heat = PH(h_ind,3);   %total heating power in W
    P_inp = (P_heat - P_rad) / 2; %mainly distributed to e and i
    P_inp = P_inp / e;  %power in eV/s

    Da = - P_inp ./ (S .* ne .* dT);

    fid = fopen(['./out/', num2str(shot), '_', num2str(time), '_Da.dat'], 'w');
    fprintf(fid, '%g\t%g\n', [r';Da']);
    fclose(fid);
    
    %##########################################################################
    % PLOT
    %##########################################################################

    figure('units', 'normalized', 'outerposition', [0, 0, 1, 1]);
    axis tight
    semilogy(r, Da, '--r', 'LineWidth', 2, 'DisplayName', 'D_a estimated');
    hold on
    semilogy(r, 1e4 .* ones(size(r)), '-k', 'LineWidth', 2, 'DisplayName', 'D_a = 10^4 cm^2 s^{-1}');
    hold off
    xlim([0, max(r)])
    ylim([1 1e8])
    legend('Location', 'best')
    xlabel('r / cm')
    ylabel('D / cm^2 s^{-1}')
    title('Estimation of D_a with total heating power')
    set(findall(gcf,'-property','FontSize'), 'FontSize', 18)

    print(['./out/', num2str(shot), '_', num2str(time), '_Da.png'], '-dpng', '-r200')
    print(['./out/', num2str(shot), '_', num2str(time), '_Da.svg'], '-dsvg')
    
    close all
    
    %##########################################################################
    % PLOT S
    %##########################################################################

    % figure('units', 'normalized', 'outerposition', [0, 0, 1, 1]);
    % axis tight
    % plot(r, S, '--r', 'LineWidth', 2, 'DisplayName', 'S from equil r q psi');
    % hold on
    % plot(r, S_est, '-k', 'LineWidth', 2, 'DisplayName', 'S = 4\pi^2Rr');
    % hold off
    % legend('Location', 'best')
    % xlabel('r / cm')
    % ylabel('S / cm^2')
    % set(findall(gcf,'-property','FontSize'), 'FontSize', 18)
end