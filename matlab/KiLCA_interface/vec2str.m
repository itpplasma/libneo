function s = vec2str(v, f)
%##############################################################
%function s = vec2str(v, f)
%##############################################################
% description:
%--------------------------------------------------------------
% converts a vector into a string: (a, b, ...)
%##############################################################
% input:
%--------------------------------------------------------------
% v     ... vector (row or column)
% f     ... number format according to MATLAB standard
%##############################################################
% output:
%--------------------------------------------------------------
% s     ... string
%##############################################################

%author:   Philipp Ulbl
%created:  05.02.2019
%modified: 21.08.2019
    
    %check for vector
    if ~isvector(v) || iscell(v)  %isvector is true also for cells
        error('v must be a vector.')
    end

    %number of elements
    num = numel(v);

    %create format string for sprintf
    form = '';
    for k = 1:num
        %add %g,%d,... (given by f) for each element in vector
        form = [form, f];
        %if not the last element add a , between elements
        if k < num
            form = [form, ','];
        end
    end
    
    %write vector to string
    s = sprintf(['(', form, ')'], v);        
end