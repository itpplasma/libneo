%##########################################################################
% script_GPEC.m
%##########################################################################
% description:
%--------------------------------------------------------------------------
% this script runs gpec with the MATLAB interface.
%##########################################################################

%author:   Philipp Ulbl
%created:  01.04.2020

studyname = 'BifurcThres';

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

for k = 1:size(run, 1)

    shot = run(k, 1);
    time = run(k, 2);

    runpath = ['/temp/ulbl_p/GPEC/', studyname, '/', num2str(shot), '_', num2str(time),'/'];
    gfile = ['/temp/ulbl_p/AUG/SHOTS/',num2str(shot),'/g',num2str(shot),'.',num2str(time),'_EQH'];
    cfile = ['/temp/ulbl_p/AUG/SHOTS/',num2str(shot),'/',num2str(shot),'.',num2str(time), '_coil.dat'];

    gpec = GPEC_interface(runpath, shot, time, studyname);
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