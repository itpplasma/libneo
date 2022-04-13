function Kisslinger(coilpath, shot, cfile, pfile)
%##########################################################################
%function Kisslinger(coilpath, shot, cfile, pfile)
%##########################################################################
% description:
%--------------------------------------------------------------------------
% runs Kisslinger.x located in coilpath for shot and cfile and saves
% field.dat in pfile.
%
% uses coil file from experiment which contains 
% Il1, Il2, ... Il8, Iu1, Iu2, ... Iu8 in a row and converts this to the
% proper format of cur_asd.dd.
%##########################################################################
% input:
%--------------------------------------------------------------------------
% coilpath ... location of kisslinger.x
% shot     ... shot number
% cfile    ... coil file from experiment to use
% pfile    ... path to field.dat to generate
%##########################################################################
            
    if(exist(cfile, 'file') ~= 2)
        error(['cfile not found in: ', cfile]);
    end
    
    mpath = pwd();
    
    %COIL CONVERTER
    curfile = [coilpath, 'cur_asd/cur_asd_', num2str(shot), '.dd'];
    coil = AUG_coil(cfile);
    coil.read(false);
    coil.export2Kisslinger(curfile, false);

    cd(coilpath)
    %create symbolic link to pfile path
    [ppath, ~, ~] = fileparts(pfile);
    ppath = [ppath, '/'];
    %create path to pfile if not exists
    system(['mkdir -p ', ppath]);
    %create symbolic link to pfile
    system(['ln -sfT ', ppath, ' ', 'MESH3D']);
    %create symbolic link to cur file
    system(['ln -sf ', curfile, ' ', 'cur_asd.dd']);
    %KISSLINGER.INP - could be set here
    
    %run kisslinger to create field.dat. ~1.5h run time
    time_start = datestr(datetime);
    disp(['Start of Kisslinger at ', time_start])
    [~, res] = system('./kisslinger_asdex.x');
    time_end = datestr(datetime);
    
    %log file
    logfile = [coilpath, 'kisslinger_', strrep(time_end, ' ', '_'), '.log'];
    fid = fopen(logfile, 'w');
    fprintf(fid, '%s\n', ['Start of Kisslinger at ', time_start]); %write start time
    fprintf(fid, '%s\n', ['Finished Kisslinger at ', time_end]); %write end time
    fprintf(fid, '%s\n', res); %write console output
    fclose(fid);
    
    %move log to location of field.dat (for debug)
    system(['mv -f ' logfile, ' ', ppath]);
    %copy cur_asd.dd to location of field.dat (for debug)
    system(['cp -f ' curfile, ' ', ppath]);
    
    disp(['Finished Kisslinger at ', time_end])
    cd(mpath)
end
