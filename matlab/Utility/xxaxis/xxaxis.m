function b = xxaxis(apos, bpos, height)
%##########################################################################
%function b = xxaxis(apos, bpos, height)
%##########################################################################
% description:
%--------------------------------------------------------------------------
% creates a 2nd x axis below the 1st axis (gca) that already exists.
% vertical positions as well as the plot height can be adjusted but are
% also set by default values.
%##########################################################################
% input:
%--------------------------------------------------------------------------
% apos      ... y position of 1st axis (default = 0.3)
% bpos      ... y position of 2nd axis (new) (default = 0.15)
% height    ... height of the plot (1st axis) (default = 0.6)
%##########################################################################
% output:
%--------------------------------------------------------------------------
% b         ... handle to 2nd axis
%##########################################################################

%author:   Philipp Ulbl
%created:  04.03.2020

    %set default values
    if(nargin < 3 || isempty(height))
        height = 0.6;
    end
    if(nargin < 2 || isempty(bpos))
        bpos = 0.15;
    end
    if(nargin < 1 || isempty(apos))
        apos = 0.3;
    end
    
    %set 1st axis to normalized units and shift position upwards
    a=gca;
    set(a, 'Units', 'normalized');
    pos = get(a, 'Position');
    set(a, 'Position', [pos(1), apos, pos(3), height])
    
    %create second axis with height = 1e-12 below 1st axis and remove
    %color
    b=axes('Position', [pos(1), bpos, pos(3), 1e-12]);
    set(b, 'Units', 'normalized');
    set(b, 'Color', 'none');
end
