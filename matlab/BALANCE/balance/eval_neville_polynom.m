function res = eval_neville_polynom(dim, xg, yg, deg, x, Dmin, Dmax)
%##############################################################
%function res = eval_neville_polynom(dim, xg, yg, deg, x, Dmin, Dmax)
%##############################################################
% description:
%--------------------------------------------------------------
% Legacy function of Ivan Ivanov for Neville polynomial
% interpolation.
%##############################################################
% input:
%--------------------------------------------------------------
% dim   ... dimension (length of xg?)
% xg    ... given x points
% yg    ... given y points
% deg   ... degree of polynomial
% x     ... new x points for evaluation (single point)
% Dmin  ... minimum degree of derivative
% Dmax  ... maximum degree of derivative
%##############################################################
% output:
%--------------------------------------------------------------
% res   ... neville polynomial and derivatives evaluated at x
%##############################################################
    
%author:   Philipp Ulbl
%created:  10.01.2020

    %if x < xg(1) | x > xg(dim),
    %    disp ('warning: eval_neville_polynom: x is outside the grid');
    %end

    ind = find_index_for_interp (deg, x, dim, xg);

    xa = xg(ind+1 : ind+deg+1);
    ya = yg(ind+1 : ind+deg+1);

    %if x < xa(1) | x > xa(deg+1),
    %    disp ('warning: eval_neville_polynom: x is outside the interpolation array');
    %end

    p = zeros (Dmax+2, deg+1, deg+1);

    for d = 0:Dmax+1
        for i = 0:deg
            for j = 0:deg
                p(d+1,i+1,j+1) = 0.0e0;
            end
        end
    end

    for d = 0:deg
        p(2,d+1,d+1) = ya(d+1);
    end

    for n = 0:Dmax %over derivative order
        for d = 1:deg %over polymomial order
            for i = 0:deg-d
                j = i + d;

                p(n+2,i+1,j+1) = n*(p(n+1,i+1,j) - p(n+1,i+2,j+1))/(xa(i+1)-xa(j+1)) + ...
                                 ((x-xa(j+1))/(xa(i+1)-xa(j+1)))*p(n+2,i+1,j) - ...
                                 ((x-xa(i+1))/(xa(i+1)-xa(j+1)))*p(n+2,i+2,j+1);
            end
        end
    end

    %store results:
    for n = Dmin:Dmax
        res(n-Dmin+1) = p(n+2,1,deg+1);
    end
end

function ind = find_index_for_interp (deg, x, dimx, xa)

    ind = search_array (x, dimx, xa);

    ind = round (min (max (0, ind-(deg-1)/2), dimx-deg-1));
end

function ind = search_array (x, dimx, xa)

    %0 <= ind <= dimx-2
    %warning: returns dim-2 for x>=xa[dim-1] - Ok for splines and interpolation
    %be careful to change something here:

    ilo = 0;
    ihi = dimx;

    while ihi > ilo+1
        i = fix((ihi+ilo)/2);
        if xa(i+1) > x
            ihi = i;
        else
            ilo = i;
        end
    end

    ind = ilo;
end