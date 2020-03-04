function out = smooth2level(y, x, level, span, method)
%##########################################################################
%function out = smooth2level(y, x, level, span, method)
%##########################################################################
% description:
%--------------------------------------------------------------------------
% smooths level of derivative of datapoints y,x with span number/percent of
% datapoints using method. If level = 0, only datapoints are smoothed, if
% level > 0, all derivatives are smoothed as well.
% function first smooths, then computes the derivative using gradient. This
% is done until level is reached. In the end, starting from the highest
% order derivative, smoothed values are integrated back to reach the
% previous derivatives and the datapoints themselves (cumtrapz). The
% integration constant is chosen to be the first value of the previous
% derivative / datapoint.
% 
% NOTICE: smoothing to arbitrary level is not a good idea.
%         you can try with Er.dat using level=1 and level=2 or with
%         Te.dat using level=2 and level=3. datapoints themselves may get
%         off the original ones if you use too high level
%##########################################################################
% input:
%--------------------------------------------------------------------------
% y         ... y - datapoints
% x         ... x - datapoints
% level     ... level of derivative to smooth (0 = function itself)
% span      ... number of datapoints for calculating smoothed values (from
%               smooth function of matlab). 
%                 or
%               percent of datapoints for calculating smoothed values
% method    ... method for smoothing (optional, default = from smooth)
%##########################################################################
% output:
%--------------------------------------------------------------------------
% out       ... cell array with dimension of level + 1 containing smoothed
%               datapoints and derivatives
%##########################################################################

%author:   Philipp Ulbl
%created:  21.01.2020

    %construct cell array with quantity and placeholder for derivatives
    out = [{y}, cell(1, level)];
    
    %compute all derivatives up to level and smooth
    for j = 0:level
        %compute derivative
        if(j > 0)
            out{j+1} = gradient(out{j}, x);
        end
        %smoothen derivative
        if(nargin < 5 || isempty(method))
            out{j+1} = smooth(x.^2, out{j+1}, span);
        else %if method specified
            out{j+1} = smooth(x.^2, out{j+1}, span, method);
        end
    end

    %integrate back up to quantity
    for j = 1:level
        ind = level + 1 - j;
        out{ind} = cumtrapz(x, out{ind+1}) + out{ind}(1);
    end
end