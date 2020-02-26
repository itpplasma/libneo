%##########################################################################
%
% This script shows how to use the efit class to read 2 G-EQDSK files in 
% and compare them to eachother.
%
%##########################################################################

%author: Philipp Ulbl
%created: 23.05.2019

%##########################################################################
% choose a file to read
%##########################################################################

fname1 = './example/g30835.3200_ed6';
fname2 = './example/g147131.02300_DIIID_KEFIT';
fname1='/proj/plasma/RMP/DATA2017/33133/3.0s/g33133.3000_ed4';
fname2='/temp/ulbl_p/AUG/SHOTS/33133/g33133.3000_EQH';
%##########################################################################
% create an instance and use it to read the file
%##########################################################################

%create objects
e1 = efit(fname1, [], []);
e2 = efit(fname2, [], []);

%read from files
e1.read();
e2.read();

%##########################################################################
% compare both files against each other
%##########################################################################

e1.compareto(e2); %output in console + figure