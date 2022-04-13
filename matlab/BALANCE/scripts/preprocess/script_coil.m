%##########################################################################
% script_coil.m
%##########################################################################
% description:
%--------------------------------------------------------------------------
% this script runs Kisslinger.x for each specified shot
%##########################################################################
            
%author:   Philipp Ulbl
%created:  17.01.2020

%add necessary paths
libBalance = '~/BALANCE/';
addpath(libBalance)
coilpath = [libBalance, 'coil/'];

%shots and times to calculate (single time should be sufficient per shot because RMPs are static)
%shot = [33120, 33133, 33353, 34213, 34548];
%time = [ 5500,  2600,  2325,  2865,  5619];

shot = [33353, 34213, 34548];
time = [ 2325,  2865,  5619];

disp(['calculation for shots: ', num2str(shot)]);

%calculation for each shot
for k = 1:numel(shot)
    try
        %get name of experimental coil file
        filehead = ['/temp/ulbl_p/AUG/SHOTS/',num2str(shot(k)),'/',num2str(shot(k)),'.',num2str(time(k))];
        cfile  = [filehead,'_coil.dat'];
        %get directory of field.dat
        pfile = ['/temp/ulbl_p/MESH3D/',num2str(shot(k)),'/field.dat'];
        %run Kisslinger
	disp(['run ', num2str(k), '/', num2str(numel(shot)), ' shot ', num2str(shot(k))])
        Kisslinger(coilpath, shot(k), cfile, pfile);
    catch ex
        disp(ex.message)
    end
end
