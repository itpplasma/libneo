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

addpath(libKiLCA)
addpath(libBalance)
addpath([libBalance, 'balance/'])

mpath = pwd();

studyname = 'DaEstimation/';
system(['mkdir -p ~/Balance_Results/', studyname]);

%Runs to make
run = [30835, 2300;...
       30835, 3000;...
       33120, 5500;...
       33120, 5635;...
       33133, 2600;...
       33133, 3000;...
       33353, 2325;...
       33353, 2670;...
       33353, 2900];

m = 4:9;
n = 2 .* ones(size(m));

for k = 1:size(run, 1)

    shot = run(k, 1);
    time = run(k, 2);

    runpath = ['/temp/ulbl_p/BALANCE_2020/', studyname, num2str(shot), '_', num2str(time),'/'];
    filename = strsplit(runpath, '/');
    filename = filename{end-1};

    %experimental input
    gfile  = ['/temp/ulbl_p/AUG/SHOTS/',num2str(shot),'/g',num2str(shot),'.',num2str(time),'_EQH'];
    filehead = ['/temp/ulbl_p/AUG/SHOTS/',num2str(shot),'/',num2str(shot),'.',num2str(time)];
    cfile  = [filehead,'_coil.dat'];

    neprof = [filehead,'_ne_PED_ulbl_rho_pol.dat'];
    Teprof = [filehead,'_Te_PED_ulbl_rho_pol.dat'];
    Tiprof = [filehead,'_Ti_PED_ulbl_rho_pol.dat'];
    vtprof = [filehead,'_vt_PED_ulbl_rho_pol.dat'];
    
    runprof = ['/temp/ulbl_p/BALANCE_2020/suppression/', num2str(shot), '_', num2str(time), '/profiles/'];
    system(['mkdir -p ', runpath]);
    system(['rsync -a ', runprof, ' ', runpath, 'profiles/']);
        
    %location of field.dat
    pfile = ['/temp/ulbl_p/MESH3D/',num2str(shot),'/field.dat']; %will be calculated if not present
    fluxdatapath = ['/temp/ulbl_p/FLUXDATA/',num2str(shot),'/',num2str(time),'/']; %will be calculated if not present

    %##########################################################################
    % RUN BALANCE WITH MATLAB CLASS
    %##########################################################################

    %BALANCE CODE
    bal = Balance(runpath, shot, time, libKiLCA);
    bal.FLAG_FORCE_PROFILES = false;
    bal.setCoil(pfile, cfile);
    bal.setEqui(gfile, fluxdatapath);
    bal.setProfiles(neprof, Teprof, Tiprof, vtprof);
    bal.setModes(m, n);
    bal.write();
    bal.run();

    %bal.loadOutput();
    %bal.basicPost()

    %##########################################################################
    % BALANCE CURRENTS
    %##########################################################################

    Brvac = zeros(size(m));
    for l=1:numel(m)
        post = bal.kil_vac.postprocessors{l};
        Brvac(l) = interp1(post.r, post.Br, post.rres);
    end
    A_ref = (bal.r_res .* bal.r_big ./ bal.n) .* Brvac;

    bal.readAntenna();
    A_realistic = zeros(size(m));
    for l=1:numel(m)
        A_realistic(l) = interp1(bal.Atheta_antenna{l}(:, 1), bal.Atheta_antenna{l}(:, 2), bal.r_res(l));
    end

    bal.rescaleI(abs(A_realistic ./ A_ref));

    %##########################################################################
    % READ GPEC CURRENTS
    %##########################################################################

    gpecpath = ['/temp/ulbl_p/GPEC/', num2str(shot), '_', num2str(time), '/'];

    qraw = ncread([gpecpath, 'gpec_profile_output_n2.nc'], 'q_rational'); %read I res in A
    Iraw = ncread([gpecpath, 'gpec_profile_output_n2.nc'], 'I_res'); %read I res in A
    Igpec = Iraw(:, 1) + 1i .* Iraw(:, 2);
    ind = ismember(abs(qraw), abs(m./n));
    Igpec = abs(Igpec(ind)') / 10; %cgs

    bal.rescaleD((Igpec./bal.I_res_resc).^2);

    %##########################################################################
    % TABLE OF CURRENTS
    %##########################################################################

    colnames = {'SpecFac', 'I_KiLCA', 'I_rescaled', 'I_GPEC', 'FormFac', 'De22', 'De22_rescaled'};
    rownames = arrayfun(@(x,y) [num2str(x), ',', num2str(y)], m, n, 'UniformOutput', false);
    tableval = {bal.specfactors', bal.I_res', bal.I_res_resc', ...
                Igpec', bal.formfactors', bal.De22_res', bal.De22_res_resc'};

    fid = fopen(['~/Balance_Results/', studyname, filename, '_CurTable.txt'], 'w');

    row = 'm,n\t';
    for j = 1:numel(colnames)
        row = [row, pad(colnames{j}, 12), '\t'];
    end
    fprintf(fid, [row, '\n']);

    for i = 1:numel(rownames)
        row = [rownames{i}, '\t'];
        for j = 1:numel(colnames)
            row = [row, sprintf('%e', tableval{j}(i)), '\t'];
        end
        fprintf(fid, [row, '\n']);
    end
    fclose(fid);
end

%##########################################################################
% PROCESS
%##########################################################################

% figure('units', 'normalized', 'outerposition', [0, 0, 1, 1]);
% axis tight
% bal.plotProfiles();
% print(['~/Balance_Results/', studyname, filename, '_Prof'], '-dpng', '-r200')
% %print(['~/Balance_Results/', studyname, filename, '_Prof'], '-dsvg')
% 
% figure('units', 'normalized', 'outerposition', [0, 0, 1, 1]);
% axis tight
% bal.plotDe22Res();
% print(['~/Balance_Results/', studyname, filename, '_De22Res'], '-dpng', '-r200')
% print(['~/Balance_Results/', studyname, filename, '_De22Res'], '-dsvg')








