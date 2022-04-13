%##########################################################################
% function save_file(data, name)
%##########################################################################
% DESC 
%--------------------------------------------------------------------------
% saves a file with specified name containing given data
%##########################################################################
% INPUT
%--------------------------------------------------------------------------
% data      ... can be matrix or cell-array of strings
% name      ... name of the file (with path)
%##########################################################################

%author:   Philipp Ulbl
%created:  04.02.2019
%modified: 21.08.2019

function save_file(data, name)

    %open file to write
    fid = fopen(name, 'w');
    
    %save cell array of strings
    if iscell(data) && isvector(data)
        fprintf(fid, '%s\n', data{:});
    %or matrix
    elseif ~iscell(data) && ismatrix(data)
        %build format string for rows
        form = '';
        for k = 1:size(data, 2)
            form = [form, '%f '];
        end
        %fprint writes row-wise into file, column format given by form
        fprintf(fid, [form, '\n'], data');
    else
        error('datatype not supported (only matrix or cell-array of strings).')
    end
    %close file
    fclose(fid);
end