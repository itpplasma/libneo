function raw = read_in(path)
%##############################################################
%function raw = read_in(path)
%##############################################################
% description:
%--------------------------------------------------------------
% reads the file specified by path and returns the raw data in
% form of a cell array which contains each row in the file as
% rows in the cell array
%##############################################################
% input:
%--------------------------------------------------------------
% path  ... path of file
%##############################################################
% output:
%--------------------------------------------------------------
% raw   ... raw data as 1d row-wise cell array
%##############################################################

%author:   Philipp Ulbl
%created:  05.02.2019
%modified: 05.02.2019
    
    %open file
    fid = fopen(path);
    %get a line (initialize tline)
    tline = fgetl(fid);
    raw = {};
    
    %continue while tline is not empty
    while ischar(tline)
        %add tline
        raw(end+1) = {tline};
        %and get new line
        tline = fgetl(fid);
    end
    raw = raw';
    
    %close file
    fclose(fid);
end