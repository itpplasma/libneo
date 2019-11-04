%##########################################################################
%
% This script shows how to use the mnDAT class to create mndat files in 
% MATLAB.
%
%##########################################################################

%author: Philipp Ulbl
%created: 04.11.2019

%##########################################################################
% create some arbitrary data
%##########################################################################

n = 2;
mmax = 3;
m = [-mmax:-1, 1:mmax];

%create flux surface label
s = linspace(0, 1, 10)';
%create q
q = 1 + s.^2 .* 6;
%create radial B
Br = s * sqrt(m);
%create mn matrix
mn = [real(Br), -imag(Br)];

%##########################################################################
% create an instance of the class and use it
%##########################################################################

%construct class with path
mndat = mnDAT('./out/', 'B', 'r', n, true);

%set parameters
mndat.set(s, q, mn);

%write kin file
mndat.write();