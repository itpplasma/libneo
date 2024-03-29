%##########################################################################
% script_balance.m
%##########################################################################
% description:
%--------------------------------------------------------------------------
% this script runs the balance code using function singlerun_balance.
% it requires data to be present in /temp/ulbl_p/ directory (see function)
%##########################################################################

%author:   Philipp Ulbl
%created:  25.03.2020

clear all

libKiLCA = '~/KiLCA_interface/';
libBalance = '~/BALANCE/balance';

addpath(genpath(libKiLCA))
addpath(genpath(libBalance))

mpath = pwd();

studyname = 'Simple_33120_5500_Test';
system(['mkdir -p ~/Balance_Results/', studyname, '/']);

%Runs to make
% run = [30835, 2300;...
%        30835, 3000;...
%        33120, 5500;...
%        33120, 5635;...
%        33133, 2600;...
%        33133, 3000;...
%        33353, 2325;...
%        33353, 2670;...
%        33353, 2900];

run = [33120, 5500];
%copy = '/temp/ulbl_p/BALANCE_2020/TimeEvol_ReprMartin33133/33133_3000/profiles/';
copy = '';

m = 4:9;
n = 2.*ones(size(m));
modes = [m' , n'];
   
for k = 1:size(run, 1)

    shot = run(k, 1);
    time = run(k, 2);

    %gpecpath = ['/temp/ulbl_p/GPEC/TimeEvol_ReprMartin33133/', num2str(shot), '_', num2str(time), '/'];
    %gpecpath = ['/temp/ulbl_p/GPEC/', studyname, '/', num2str(shot), '_', num2str(time), '/'];
    gpecpath = '/temp/ulbl_p/GPEC/test1/';
    
    bal = singlerun_balance(studyname, shot, time, false, modes, copy, gpecpath);
end
