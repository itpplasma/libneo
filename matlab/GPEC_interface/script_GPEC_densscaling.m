%##########################################################################
% script_GPEC_densscaling.m
%##########################################################################
% description:
%--------------------------------------------------------------------------
% this script runs gpec with the MATLAB interface. density scaling runs
%##########################################################################

%author:   Philipp Ulbl
%created:  08.04.2020

addpath('../EFIT/')

studyname = 'Density_Study3';

%Runs to make
shot = 33353;
time = 2900;

runpath = ['/temp/ulbl_p/GPEC/', studyname, '/', num2str(shot), '_', num2str(time),'/'];
gfile = ['/temp/ulbl_p/AUG/SHOTS/',num2str(shot),'/g',num2str(shot),'.',num2str(time),'_EQH'];
cfile = ['/temp/ulbl_p/AUG/SHOTS/',num2str(shot),'/',num2str(shot),'.',num2str(time), '_coil.dat'];

nfac = [0.1, 0.325, 0.55, 0.775, 1, 2, 4, 6, 8, 10];

for k = 1:numel(nfac)

    e = efit(gfile, [], []);
    e.read();
    e.pres = e.pres .* nfac(k);
    e.pprime = e.pprime .* nfac(k);
    e.write(['n', num2str(nfac(k))]);
    gfile = e.fname;
    
    gpec = GPEC_interface([runpath, num2str(nfac(k)), '/'], shot, time, studyname);
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

col = [166,206,227;...
31,120,180;...
178,223,138;...
51,160,44;...
251,154,153;...
227,26,28;...
253,191,111;...
255,127,0;...
202,178,214;...
106,61,154];
col = col./255;

figure('units', 'normalized', 'outerposition', [0,0,1,1])
for k = 1:numel(nfac)
    subplot(2,1,1)
    fname = [runpath, num2str(nfac(k)), '/gpec_profile_output_n2.nc'];
    q = ncread(fname, 'q_rational');
    I = ncread(fname, 'I_res');
    I = abs(I(:, 1) + 1i .* I(:, 2));
    
    p=semilogy(q*2, I, '--o', 'LineWidth', 2, 'DisplayName', ['n = ', num2str(nfac(k)), 'n_0'], 'Color', col(k, :));
    set(p, 'MarkerFaceColor', get(p, 'Color'));
    hold on
end
hold off
legend()
xlabel('m')
ylabel('|I^{res}_{mn}| / A')
for k = 1:(numel(nfac)-3)
    subplot(2,1,2)
    fname = [runpath, num2str(nfac(k)), '/gpec_profile_output_n2.nc'];
    q = ncread(fname, 'q_rational');
    I = ncread(fname, 'I_res');
    I = abs(I(:, 1) + 1i .* I(:, 2));
    
    p=semilogy(q*2, I, '--o', 'LineWidth', 2, 'DisplayName', ['n = ', num2str(nfac(k)), 'n_0'], 'Color', col(k, :));
    set(p, 'MarkerFaceColor', get(p, 'Color'));
    hold on
end
hold off
legend()
xlabel('m')
ylabel('|I^{res}_{mn}| / A')
set(findall(gcf,'-property','FontSize'), 'FontSize', 16)
print('~/Meetings/Meeting_2020_04_09/GPEC_currents.png', '-dpng', '-r200')