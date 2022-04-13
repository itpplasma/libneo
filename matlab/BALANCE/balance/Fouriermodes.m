function Fouriermodes(fourierpath, fluxdatapath, fdb0)
%##########################################################################
%function Fouriermodes(fourierpath, fluxdatapath, fdb0)
%##########################################################################
% description:
%--------------------------------------------------------------------------
% runs fouriermodes.x located in fourierpath for already written
% field_divB0.inp file (MATLAB object fdb0) and saves output to
% fluxdatapath.
%##########################################################################
% input:
%--------------------------------------------------------------------------
% fourierpath  ... location of fouriermodes.x
% fluxdatapath ... location to save data
% fdb0         ... field_divB0.m object (already written)
%##########################################################################
            
    mpath = pwd();

    cd(fourierpath)
    %FOURIERMODES.INP - could be set here
    time_start = datestr(datetime);
    disp(['Start of Fouriermodes at ', time_start])
    [stat, res] = system('./fouriermodes.x');
    time_end = datestr(datetime);

    %log file
    logname = [fourierpath, 'fouriermodes_', strrep(time_end, ' ', '_'), '.log'];
    fid = fopen(logname, 'w');
    fprintf(fid, '%s\n', ['Start of Fouriermodes at ', time_start]); %write start time
    fprintf(fid, '%s\n', ['Finished Fouriermodes at ', time_end]); %write end time
    fprintf(fid, '%s\n', res); %write console output
    fclose(fid);

    if(stat ~= 0)
        error('error in fouriermodes, see log file.')
    end
    
    %move output to fluxdatapath
    system(['mkdir -p ', fluxdatapath]);
    outfiles = {'amn.dat', 'btor_rbig.dat', 'equil_r_q_psi.dat', ...
                'axis.dat', 'box_size.dat', 'separ.dat', 'phinorm_arr.dat', ...
                'thetabooz.dat', 'theta_of_theta_qt_flabel.dat', logname};
    for k=1:numel(outfiles)
        %skip amn.dat if run was only equi without coils
%         if(fdb0.ipert == 0 && strcmp(outfiles{k}, 'amn.dat'))
%             system(['rm ' outfiles{k}]);
%             continue;
%         end
    	system(['mv -f ' outfiles{k}, ' ', fluxdatapath]);
    end
    
    %copy input files to fluxdatapath (for debug)
    system(['cp -f ' 'field_divB0.inp', ' ', fluxdatapath]);
    system(['cp -f ' 'fouriermodes.inp', ' ', fluxdatapath]);

    disp(['Finished Fouriermodes at ', time_end])
    
    cd(mpath)
end
