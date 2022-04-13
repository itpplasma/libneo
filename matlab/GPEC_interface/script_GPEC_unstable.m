%##########################################################################
% script_GPEC.m
%##########################################################################
% description:
%--------------------------------------------------------------------------
% this script runs gpec with the MATLAB interface.
%##########################################################################

%author:   Philipp Ulbl
%created:  18.05.2020

studyname = 'UnstableRuns';

%Runs to make
run = [30835, 3000];

for k = 1:size(run, 1)

    shot = run(k, 1);
    time = run(k, 2);

    runpath = ['/temp/ulbl_p/GPEC/', studyname, '/', num2str(shot), '_', num2str(time),'/'];
    gfile = ['/temp/ulbl_p/AUG/SHOTS/',num2str(shot),'/g',num2str(shot),'.',num2str(time),'_EQH'];
    cfile = ['/temp/ulbl_p/AUG/SHOTS/',num2str(shot),'/',num2str(shot),'.',num2str(time), '_coil.dat'];

    gpec = GPEC_interface(runpath, shot, time, studyname);
    gpec.setEqui(gfile);
    gpec.setCoil(cfile);
    gpec.loadDGV();
    gpec.equil.EQUIL_CONTROL.psihigh = 0.993; %unstable value for this shot/time = 0.992
    %gpec.dcon.DCON_CONTROL.peak_flag = true;
    gpec.dcon.DCON_CONTROL.sas_flag = false;
    gpec.dcon.DCON_CONTROL.psiedge = 0.988;
    gpec.write();
    
    try
        gpec.run();
    catch
        warning('error in GPEC occurred.')
        gpec.run();
    end
        
end



