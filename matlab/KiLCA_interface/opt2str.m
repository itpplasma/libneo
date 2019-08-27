function s = opt2str(o)
%##############################################################
%function s = opt2str(o)
%##############################################################
% description:
%--------------------------------------------------------------
% converts an option into a string. option is meant in the 
% context of KiLCA input parameters and the conversion is done
% for the format used in KiLCA.
%##############################################################
% input:
%--------------------------------------------------------------
% o     ... 'option': char, scalar, vector
%##############################################################
% output:
%--------------------------------------------------------------
% s     ... string
%##############################################################

%author:   Philipp Ulbl
%created:  05.02.2019
%modified: 21.08.2019
    
    %check for datatype and convert to string
    if ischar(o)
        s = o;
    elseif isscalar(o)
        s = num2str(o, '%g');
    elseif isvector(o) && ~iscell(o) %isvector is true also for cells
        s = vec2str(o, '%g');
%     elseif iscell(o) %maybe change??? not used by now...
%         s = o;
    else
        error('type of option not supported (only char, scalar, vector).')
    end
end