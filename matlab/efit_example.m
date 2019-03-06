%##########################################################################
%
% This script shows how to use the efit class for MATLAB.
%
%##########################################################################

%author: Philipp Ulbl
%created: 15.02.2019
%modified: 05.03.2019

%##########################################################################
%generate data for testing the class:
%##########################################################################

npoints = 1e2;

%parameters
rmin = 1e-3;
r0 = 0.67;
q0 = 8;
T0 = 8e3;
n0 = 1e13;

%r profiles and s = f(r)
r = linspace(rmin, r0, npoints);
theta = linspace(0, 2*pi, npoints);
s = (r./r0).^2;

%test profiles as functions of s
q = -(1 + q0 .* s); %safety factor
T = T0 .* (1 - s); %temperature
n = n0 .* (1 - s); %num density
p = n .* T; %pressure

%calculate some psi profile
b0z_a = -1.752e4.*ones(size(r));
psi_a = (b0z_a.*((r0^2)/(2*q0))) .* log(1+q0.*(r./r0).^2);

%convert psi-profile to circular grid
[RR, THTH] = meshgrid(r, theta);
RM = RR .* cos(THTH);
ZM = RR .* sin(THTH);
PSI = psi .* sqrt((RM*5).^2+ZM.^2);

%convert circular psi to RZ grid using interpolation
R_efit = [-fliplr(r), r(2:end)];
Z_efit = [-fliplr(r), r(2:end)];
[RM_efit, ZM_efit] = meshgrid(R_efit, Z_efit);
PSI_efit = griddata(RM, ZM, PSI, RM_efit, ZM_efit);
PSI_efit(isnan(PSI_efit)) = 0;

%##########################################################################
% create an instance and use it to write the file
%##########################################################################

%recalc q and p for EFIT grid (from r0->r profile to -R->R)
Q = [fliplr(q), q(2:end)];
P = [fliplr(n), n(2:end)] .* [fliplr(Te), Te(2:end)];

%create instance of efit class in current folder (pwd)
e = efit([pwd(), '/'], 2*npoints-1, 2*npoints-1);
e.CASE{5} = 'TESTDATA';
e.CASE{6} = 'ONLYPSI';

%fill efit class with given data
e.rdim = r0;
e.zdim = r0;
e.rleft = rmin;

e.zmid = 0;
e.rmaxis = 0;
e.zmaxis = 0;

e.pres = P;
e.qpsi = Q;

e.psirz = PSI_efit;

%some arbitrary values to test the last data block
e.nbbbs = 5;
e.rbbbs = linspace(1, 5, 5);
e.zbbbs = linspace(-1, -5, 5);

%fill empty properties with zeros
e.fillempty();

%write to file
e.write();