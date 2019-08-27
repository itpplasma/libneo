function c = struct2plaincell(s)
%##############################################################
%function c = struct2plaincell(s)
%##############################################################
% description:
%--------------------------------------------------------------
% converts a struct into a plain cell array (recursive).
% plain cell array: no subcells or substruct, only scalars, 
% strings or arrays.
%##############################################################
% input:
%--------------------------------------------------------------
% s     ... struct
%##############################################################
% output:
%--------------------------------------------------------------
% c     ... plain cell array
%##############################################################

%author:   Philipp Ulbl
%created:  06.02.2019
%modified: 06.02.2019
    
    c = {};
    %convert struct to cell (may have nested structs by now)
    s = struct2cell(s);
    
    %iterate all elements of s and check the datatype
    for k = 1:numel(s)
        %if not struct write element into cell array
        if ~isstruct(s{k})
            %this writes char
            if ischar(s{k})
                c{end+1} = s{k};
            %this writes scalars, vectors and arrays
            else
                %arrays will be written row-wise
                for l = 1:size(s{k}, 1)
                    c{end+1} = s{k}(l, :);
                end
            end
        %if sub-element is struct, recursive call of this method
        else
            %convert sub-struct to plain cell
            s2 = struct2plaincell(s{k});
            %add all elements of sub-struct to this cell
            for l = 1:numel(s2)
                c{end+1} = s2{l};
            end
        end
    end
    
    %finally transpose for right dimensions
    c = c';
end