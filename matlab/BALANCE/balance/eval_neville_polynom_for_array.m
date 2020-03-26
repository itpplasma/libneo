function res = eval_neville_polynom_for_array (dim, xg, yg, deg, x, Dmin, Dmax)
%##############################################################
%function res = eval_neville_polynom(dim, xg, yg, deg, x, Dmin, Dmax)
%##############################################################
% description:
%--------------------------------------------------------------
% Legacy function of Ivan Ivanov for Neville polynomial
% interpolation for arrays. Calls eval_neville_polynom.
%##############################################################
% input:
%--------------------------------------------------------------
% dim   ... dimension (length of xg?)
% xg    ... given x points
% yg    ... given y points
% deg   ... degree of polynomial
% x     ... new x points for evaluation (array)
% Dmin  ... minimum degree of derivative
% Dmax  ... maximum degree of derivative
%##############################################################
% output:
%--------------------------------------------------------------
% res   ... neville polynomial and derivatives evaluated at x
%##############################################################
    
%author:   Philipp Ulbl
%created:  10.01.2020

    res = zeros(length(x), Dmax-Dmin+1);

    for k = 1:length(x)
        res(k,:) = eval_neville_polynom (dim, xg, yg, deg, x(k), Dmin, Dmax);
    end
end