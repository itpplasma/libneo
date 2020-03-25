function ds = splineder(s)
%##########################################################################
%function ds = splineder(s)
%##########################################################################
% description:
%--------------------------------------------------------------------------
% Returns derivative of spline object by shift of coefficients. See
% MATLAB Answers:
% https://de.mathworks.com/matlabcentral/answers/95194-how-do-i-find-the-derivative-of-a-spline-curve-in-matlab-7-9-r2009b
%##########################################################################
% input:
%--------------------------------------------------------------------------
% s         ... spline object
%##########################################################################
% output:
%--------------------------------------------------------------------------
% ds        ... spline object representing derivative of s
%##########################################################################

    %derivative
    [breaks,coefs,l,m,d] = unmkpp(s);
    ds = mkpp(breaks,repmat(m-1:-1:1,d*l,1).*coefs(:,1:m-1),d);
end
