%##########################################################################
% script_GPEC.m
%##########################################################################
% description:
%--------------------------------------------------------------------------
% this script runs gpec with the MATLAB interface.
%##########################################################################

%author:   Philipp Ulbl
%created:  01.04.2020

studyname = 'test';

shot = 33120;
time = 5500;

runpath = ['/temp/ulbl_p/GPEC/', studyname, '/', num2str(shot), '_', num2str(time),'/'];
gfile = ['/temp/ulbl_p/AUG/SHOTS/',num2str(shot),'/g',num2str(shot),'.',num2str(time),'_EQH'];
cfile = ['/temp/ulbl_p/AUG/SHOTS/',num2str(shot),'/',num2str(shot),'.',num2str(time), '_coil.dat'];

gpec = GPEC_interface(runpath, shot, time, studyname);
gpec.setEqui(gfile);
gpec.setCoil(cfile);
gpec.write();
gpec.run();