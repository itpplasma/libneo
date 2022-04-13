%##########################################################################
%
% This script shows how to use the kin class to create kin files in 
% MATLAB.
%
%##########################################################################

%author: Philipp Ulbl
%created: 16.05.2019

%##########################################################################
% create some arbitrary profiles
%##########################################################################

%create flux surface label
s = linspace(0, 1, 100);
%densities
ni = (1-s.^2) .* 0.9e10; %ion density in m^-3 (declining quadratically)
ne = (1-s.^2) .* 1.1e10; %electron density in m^-3
ti = (1-s.^2) .* 4000;   %ion temp in eV
te = (1-s.^2) .* 8000;   %electron temp in eV
we = zeros(size(s));     %electric rotational frequency

%##########################################################################
% create an instance of the class and use it
%##########################################################################

k = kin('./out/example.kin'); %can also be used without extension .kin

k.set(s, ni, ne, ti, te, we); %set all at once
%k.s = s; %can also be used for the properties individually

%check profiles with plots
% figure
% k.plotprofile('ni'); %single profile plot
figure
k.plotsummary(); %plot of all profiles

%write kin file
k.write();