function abort = collect_kprof(kpath, fname)
%##############################################################
%function abort = collect_kprof(kpath, fname)
%##############################################################
% description:
%--------------------------------------------------------------
% Uses methods of plot_k.m by Martin Heyn to read NEO-2 output,
% extract k profile and write it to a file.
%##############################################################
% input:
%--------------------------------------------------------------
% kpath ... path of NEO-2 execution
% fname ... name of output file
%##############################################################
% output:
%--------------------------------------------------------------
% abort ... number of aborted jobs
%##############################################################   

%author:   Philipp Ulbl
%created:  20.01.2020

    %library for neo2
    addpath('/proj/plasma/Neo2/Interface/Matlab/')

    mpath = pwd();
    cd(kpath);
    
    content = dir(kpath);
    %ignore . and .. in folder
    content = content(3:end);
    
    %initialize k and boozer s
    k = [];
    s = [];
    %number of aborted jobs
    abort = 0;
    
    %iterate each element in content
    for j=1:numel(content)

        %get name of element
        dirname = content(j).name;
        
        %check if is a dir but not template
        if(content(j).isdir == true && strcmp(dirname, 'TEMPLATE_DIR') == false)
           
            try
                %extract needed files
                fulltransp  = h52struct([dirname, '/fulltransp.h5']);
                neo2config  = h52struct([dirname, '/neo2_config.h5']);

                %extract boozer s
                s(end+1) = neo2config.settings.boozer_s;
                %extract k
                try
                    k(end+1) = fulltransp.k_cof;
                catch
                    ntvout = load([dirname, '/ntv_out.dat']);
                    k(end+1) = ntvout(7);
                end
            catch
                %if file could not be read, count job as aborted
                abort = abort + 1;
            end
        end
    end
    
    %sort for ascending s
    [s, ind] = sort(s);
    k = k(ind);
    
    %load M out of surfaces
    M = load('surfaces.dat');
    %get r effective
    r = interp1(M(:,1), M(:,2), s);
    
    %export r, k into file
    export = [r', k'];
    save(fname, 'export', '-ascii');
    
    cd(mpath);
end
