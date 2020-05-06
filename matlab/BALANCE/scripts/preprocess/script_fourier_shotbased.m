%##########################################################################
% script_fourier.m
%##########################################################################
% description:
%--------------------------------------------------------------------------
% this script runs fouriermodes.x for each specified shot
%##########################################################################
            
%author:   Philipp Ulbl
%created:  17.01.2020

%add necessary paths
libBalance = '~/BALANCE/';
addpath([libBalance, 'balance/']);
fourierpath = [libBalance, 'fourier/'];

%shots and times to calculate
%shot = [33120, 33120, 33133, 33133, 33353, 33353, 33353, 34213, 34548];
%time = [ 5500,  5635,  2600,  3000,  2325,  2670,  2900,  2865,  5619];

shot = 33353;
time =  2900;

disp('calculating shots and times: ')
disp(num2str([shot',time']))

%calculation for each shot
for k = 1:numel(shot)
    try
        gfile  = ['/temp/ulbl_p/AUG/SHOTS/',num2str(shot(k)),'/g',num2str(shot(k)),'.',num2str(time(k)),'_EQH'];
        pfile = ['/temp/ulbl_p/MESH3D/',num2str(shot(k)),'/field.dat'];
        convexfile = '/proj/plasma/RMP/DATA/ASDEX/convexwall.dat';
        
        fdb0 = field_divB0(gfile, pfile, convexfile, '');
        fdb0.ipert = 1;
        fdb0.write([libBalance, 'blueprints/'], fourierpath);
        
        fluxdatapath = ['/temp/ulbl_p/FLUXDATA/',num2str(shot(k)),'/',num2str(time(k)),'/'];
        system(['mkdir -p ', fluxdatapath]);	
        
        disp(['run ', num2str(k), '/', num2str(numel(shot)), ' : shot ', num2str(shot(k)), ' ', num2str(time(k))])

        Fouriermodes(fourierpath, fluxdatapath, fdb0);
    catch ex
        disp(ex.message)
    end
end
