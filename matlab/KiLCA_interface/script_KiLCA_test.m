%##########################################################################
% THIS SCRIPT IS HERE TO QUICKLY TEST FUNCTIONS OF THE KiLCA_interface.
% -------------------------------------------------------------------------
% DO NOT USE AS A TUTORIAL, THERE IS A SPECIFIC FILE FOR THAT CAUSE!
%##########################################################################

%author:   Philipp Ulbl
%created:  06.02.2019
%modified: 25.03.2019

PLOTPROF = false;

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

path = [pwd(), '/test_newversion/'];
type = 'flre';

kil = KiLCA_interface(path, type);

kil.set_ASDEX(1);
m = 3; n = 2;
kil.modes.set(m, n);

kil.write();
kil.run();

%##########################################################################
% PLOT BACKGROUND
%##########################################################################

kil.backgrounddata.plotB();
print('./test.old/out/background_B.svg', '-dsvg')

%##########################################################################
% PLOT LINEARDATA
%##########################################################################

kil.lineardata{1}.plotB();
print('./test.old/out/linear_B.svg', '-dsvg')