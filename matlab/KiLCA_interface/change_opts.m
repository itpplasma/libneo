function raw = change_opts(raw, ind, dat)
%##############################################################
%function raw = change_opts(raw, ind, dat)
%##############################################################
% description:
%--------------------------------------------------------------
% changes in the cell array raw at the row specified by ind
% the old numerical value to a new one given by dat.
% *works also with ind and dat as arrays of same size
%##############################################################
% input:
%--------------------------------------------------------------
% raw   ... raw data
% ind   ... index to be modified
% dat   ... new data at ind
%##############################################################
% output:
%--------------------------------------------------------------
% raw   ... raw data with option modified
%##############################################################

%author:   Philipp Ulbl
%created:  05.02.2019
%modified: 21.08.2019
    
    %check for same size of dat and ind
    if size(ind) ~= size(dat)
        error('ind and dat must have the same size.') 
    end

    %go through all options to change
    for k = 1:numel(ind)
        %split the row in raw by # into option and column
        c = strsplit(raw{ind(k)}, '#');
        %change the option in the first part
        c{1} = opt2str(dat{k});
        %add whitespace between number and comment for readability
        c{2} = ['          #', c{2}];
        %save changes to row
        raw{ind(k)} = [c{:}];
    end
end