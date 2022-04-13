%##########################################################################
% THIS SCRIPT IS HERE TO QUICKLY TEST FUNCTIONS OF THE KiLCA_interface.
% -------------------------------------------------------------------------
% DO NOT USE AS A TUTORIAL, THERE IS A SPECIFIC FILE FOR THAT CAUSE!
%##########################################################################

%author:   Philipp Ulbl
%created:  06.02.2019
%modified: 25.03.2019

PLOTPROF = false;

%add class mndat to test export
addpath('../mnDAT/');

%##########################################################################
% PLOT PROFILES
%##########################################################################

names = {'q', 'Te', 'Ti', 'n', 'Vth', 'Vz', 'Er'};

if PLOTPROF == true
    figure('units', 'normalized', 'outerposition', [0, 0, 1, 1])
    
    for k = 1:numel(names)
        
        fname = ['./profiles/', names{k}, '.dat'];
        raw = importdata(fname);
        
        subplot(2,4,k)
        plot(raw(:, 1), raw(:, 2), 'LineWidth', 2)
        title(names{k})
        xlabel('r / cm')
    end
    
    print('./out/ASDEX_profiles.svg', '-dsvg')
end

%##########################################################################
% SETUP KiLCA & RUN
%##########################################################################

path = '~/KiLCA/test2/interface/';
type = 'flre';

kil = KiLCA_interface(path, type);
kil.PROF_PATH = '~/TCFP/profiles/';

m = -22; n = 2 .* ones(size(m));
kil.set_ASDEX(numel(m));
kil.modes.set(m, n);

kil.background.vgalsys = 1e9;

%kil.antenna.I0 = 1e13;
%kil.antenna.flab(1) = 1;
%
% kil.background.Btor = 23176.46;
% kil.background.flag_recalc = -1;
% kil.background.vgalsys = 1e9;
% 
% kil.zones{1}.r1 = 1e-1;
% kil.zones{1}.modelvers = 1;
% 
% kil.zones{2}.vacuum.relacc = 1e-12;
% kil.zones{2}.vacuum.absacc = 1e-12;
% kil.zones{2}.vacuum.sparse_relacc = 1e-10;
% kil.zones{2}.vacuum.sparse_absacc = 1e-10;
% 
% kil.zones{3}.vacuum.relacc = 1e-12;
% kil.zones{3}.vacuum.absacc = 1e-12;
% kil.zones{3}.vacuum.sparse_relacc = 1e-10;
% kil.zones{3}.vacuum.sparse_absacc = 1e-10;
% kil.zones{3}.vacuum.rgrid_maxdim = 100;

kil.write();
kil.run();

%##########################################################################
% PLOT BACKGROUND
%##########################################################################

kil.backgrounddata.plotB();

%##########################################################################
% PLOT LINEARDATA
%##########################################################################

if strcmp(type, 'flre')
    %kil.lineardata{1}.plotB();
    
    kil.post(1);
    %kil.postprocessors{1}.plotJcyl();
    
    %kil.Export2mnDAT('Br');
end