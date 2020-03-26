%##########################################################################
% script_balance.m
%##########################################################################
% description:
%--------------------------------------------------------------------------
% 
%##########################################################################

%author:   Philipp Ulbl
%created:  07.01.2020

libKiLCA = '~/KiLCA_interface/';
libBalance = '~/BALANCE/';

addpath(genpath(libKiLCA))
addpath(genpath(libBalance))

mpath = pwd();

studyname = 'test';
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
run = [33353, 2900];

m = 4:9;
n = 2 .* ones(size(m));

for k = 1:size(run, 1)

    shot = run(k, 1);
    time = run(k, 2);

    runpath = ['/temp/ulbl_p/BALANCE_2020/', studyname, '/', num2str(shot), '_', num2str(time),'/'];
    
    %copy profiles from suppression
    runprof = ['/temp/ulbl_p/BALANCE_2020/suppression/', num2str(shot), '_', num2str(time), '/profiles/'];
    system(['mkdir -p ', runpath]);
    system(['rsync -a ', runprof, ' ', runpath, 'profiles/']);
        
    singlerun_balance(studyname, shot, time);
end