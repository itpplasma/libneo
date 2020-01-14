%##########################################################################
%
% This script shows how to use the efit class to read G-EQDSK files in 
% MATLAB.
%
%##########################################################################

%author: Philipp Ulbl
%created: 14.05.2019

%##########################################################################
% choose a file to read
%##########################################################################

fname = './example/g33120.5500';
%fname = './example/g30835.3200_ed6';
%fname = './example/g000001.0001_TCFP';
%fname = './example/g147131.02300_DIIID_KEFIT';

%##########################################################################
% create an instance and use it to read the file
%##########################################################################

e = efit(fname, [], []);

%read from file
e.read();

%plot a single profile (1d)
%figure
%e.plot1d('q')

%plot 2d of psi-contours
%figure
%e.plot2d();

%plot OMFIT-style efit summary
e.plotsummary();