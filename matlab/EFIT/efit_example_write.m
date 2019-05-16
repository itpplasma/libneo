%##########################################################################
%
% This script shows how to use the efit class to write a G-EQDSK file in 
% MATLAB.
%
%##########################################################################

%author: Philipp Ulbl
%created: 15.02.2019

%##########################################################################
%generate data for testing the class:
%##########################################################################

%function at the end of file returns some r-vector with profiles q, p, f
%psi is given as a matrix
[r, psi, q, p, f] = efit_example_gentestdata();

%get the derivatives of p and f
pp = gradient(r, p);
ff = f .* gradient(r, f);

%##########################################################################
% create an instance and use it to write the file
%##########################################################################

n = numel(r);

%create instance of efit class in subfolder "out"
e = efit('./out/', n, n);
%add the name of the file in the header line
e.buildcase('TEST', [], [], []);

%dimension of r and z -> 2-times the 1-sided vector r
e.rdim = 2 * r(end);
e.zdim = 2 * r(end);
e.rleft = r(1); %minimum r-value

%the magnetic/geometric axis will be shifted by R=r(end) to the right.
%this changes the coordinate system from -R,R to 0,2R.

e.zmid = 0; %location of z=0 axis
e.rmaxis = r(end); %r-locatin of magnetic axis
e.zmaxis = 0; %z-location of magnetic axis

e.rcentr = r(end); %r-location of centrum (geometric axis)
e.bcentr = -1.7563e+04; %B at rcentr: some value from ASDEX settings. in G

ind = ceil(size(psi, 1)/2); %index of the middle row/col of psi
%(the r-profile of psi is located in the center row of the matrix)
e.simag = psi(ind, ind); %psi at magnetic axis
e.sibry = psi(ind, end); %psi at the border
e.current = 4.5000e+12; %current: some value from ASDEX settings. in A

%now write the profiles f, p, ff, pp, q and the matrix psi
e.fpol = f;
e.pres = p; 
e.ffprim = ff;
e.pprime = pp;
e.qpsi = q;
e.psirz = psi;

%fill empty variables with zeros
e.fillempty();

%write to file
e.write();

%##########################################################################
% function to generate testdata
%##########################################################################

function [r, psi, q, p, f] = efit_example_gentestdata()

    npoints = 129;

    %parameters
    rmin = 1e-3; %m
    r0 = 0.67; %m
    q0 = 8;
    T0 = 8e3; %kev
    n0 = 1e13; %m^-3

    %r profiles and s = f(r)
    r = linspace(rmin, r0, npoints);
    s = (r./r0).^2;

    %calculate some psi profile
    b0z_a = -1.752e4.*ones(size(r)); %G
    psi_a = (b0z_a.*((r0^2)/(2*q0))) .* log(1+q0.*(r./r0).^2);

    %test profiles as functions of s
    q = -(1 + q0 .* s); %safety factor
    T = T0 .* (1 - s); %temperature
    n = n0 .* (1 - s); %num density
    p = n .* T; %pressure
    f = r .* b0z_a;
    
    theta = linspace(0, 2*pi, npoints);

    %r-theta-circular-grid and psi based to this
    [RR, THTH] = meshgrid(r, theta);
    RM = RR .* cos(THTH);
    ZM = RR .* sin(THTH);
    PSI = psi_a .* sqrt(RM.^2+ZM.^2);

    %for efit: equi-distant RZ grid
    re = [-fliplr(r(1:2:end)), r(2:2:end)];
    [R, Z] = meshgrid(re, re);
    psi = griddata(RM, ZM, PSI, R, Z);

end