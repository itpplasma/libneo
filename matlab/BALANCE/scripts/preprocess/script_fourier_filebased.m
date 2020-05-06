%##########################################################################
% script_fourier_filebased.m
%##########################################################################
% description:
%--------------------------------------------------------------------------
% this script runs fouriermodes.x for each specified gfile
%##########################################################################
            
%author:   Philipp Ulbl
%created:  06.05.2020

%add necessary paths
libBalance = '~/BALANCE/';
addpath([libBalance, 'balance/']);
fourierpath = [libBalance, 'fourier/'];

gfiles = {'/temp/ulbl_p/AUG/WLS/33353/g33353.2900_AUGD_IDE_ed0',...
          '/temp/ulbl_p/AUG/WLS/33353/g33353.2900_MICDU_EQB_ed0',...
          '/temp/ulbl_p/AUG/WLS/33353/g33353.2900_MICDU_EQB_ed1',...
          '/temp/ulbl_p/AUG/WLS/33353/g33353.2900_MICDU_EQB_ed2',...
          '/temp/ulbl_p/AUG/WLS/33353/g33353.2900_MICDU_EQB_ed3'};

disp('calculating files: ')
cellfun(@disp, gfiles)

%calculation for each shot
for k = 1:numel(gfiles)
    
    try
        gfile  = gfiles{k};
        
        %extract shot and time from name -> assumes standard naming
        %convention
        [~, fshot, ftime] = fileparts(gfile);
        shot = str2double(strrep(fshot, 'g', ''));
        ftime = strsplit(ftime, '_');
        ftype = [ftime{2:end}];
        time = str2double(strrep(ftime{1}, '.', ''));
        
        pfile = ['/temp/ulbl_p/MESH3D/',num2str(shot(k)),'/field.dat'];
        convexfile = '/proj/plasma/RMP/DATA/ASDEX/convexwall.dat';
        
        fdb0 = field_divB0(gfile, pfile, convexfile, '');
        fdb0.ipert = 1;
        fdb0.write([libBalance, 'blueprints/'], fourierpath);
        
        fluxdatapath = ['/temp/ulbl_p/FLUXDATA/',num2str(shot(k)),'_',ftype,'/',num2str(time(k)),'/'];
        system(['mkdir -p ', fluxdatapath]);	
        
        disp(['run ', num2str(k), '/', num2str(numel(gfiles)), ' : gfile ', gfile])

        Fouriermodes(fourierpath, fluxdatapath, fdb0);
    catch ex
        disp(ex.message)
    end
end
