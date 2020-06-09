%##########################################################################
% script_parameterstudy.m
%##########################################################################
% description:
%--------------------------------------------------------------------------
% this script runs the balance code over a range of rescaled Te and n.
% legacy script
%##########################################################################
     
%author:   Philipp Ulbl
%created:  07.01.2020

libKiLCA = '~/KiLCA_interface/';
libBalance = '~/BALANCE/';

addpath(libKiLCA)
addpath(libBalance)
addpath([libBalance, 'balance/'])

mpath = pwd();

studyname = 'DensTemp_Study/';
system(['mkdir -p ~/Balance_Results/', studyname]);

%##########################################################################
% SHOT PARAMETERS
%##################/temp/ulbl_p/AUG/Suttrop/Shots/33353/2900/33353.2900_Te_PED.dat########################################################

shot = 33353;
time = 2900;

m = 4:9;
n = 2 .* ones(size(m));

%location of field.dat
pfile = ['/temp/ulbl_p/MESH3D/',num2str(shot),'/field.dat']; %will be calculated if not present
fluxdatapath = ['/temp/ulbl_p/FLUXDATA/',num2str(shot),'/',num2str(time),'/']; %will be calculated if not present

%##########################################################################
% RUN BALANCE WITH MATLAB CLASS
%##########################################################################

nfac = [linspace(0.1,1,4), linspace(1,3,9)];
nfac = unique(nfac);
Tfac = [linspace(0.5,1,5), linspace(1,2,21)];
Tfac = unique(Tfac);

for o = 1:numel(nfac)

    for p = 1:numel(Tfac)

        runpath = ['/temp/ulbl_p/BALANCE_2020/', studyname, num2str(shot), '_', num2str(time), ...
            '_n', num2str(nfac(o)), '_T', num2str(Tfac(p)),'/'];
        filename = strsplit(runpath, '/');
        filename = filename{end-1};
	
        %this is if runs were somehow disrupted to skip already calculated
        %ones
        
%         checkfile = ['out/f_', num2str(m(end)), '_', num2str(n(end)), '/fort.5000'];
%         if(exist([runpath, checkfile], 'file') == 2)
%             continue;
%         end
%         if(exist(['~/Balance_Results/DensTemp_Study/', num2str(shot), '_', num2str(time), '_n', num2str(nfac(o)), ...
%             '_T', num2str(Tfac(p)), '_CurTable.txt'], 'file'))
%             continue;
%         end

        %get gfile from proj plasma RMP
        gfile  = ['/temp/ulbl_p/AUG/SHOTS/',num2str(shot),'/g',num2str(shot),'.',num2str(time),'_EQH'];

        %get profiles from martins temp
        runprof = ['/temp/ulbl_p/BALANCE_2020/suppression/', num2str(shot), '_', num2str(time), '/profiles/'];

        %system(['mkdir -p ', runpath]);
        %system(['rsync -a ', runprof, ' ', runpath, 'profiles/']);

        %change density profile
        ne = load([runpath, 'profiles/n.dat']);
        ne(:, 2) = ne(:, 2) .* nfac(o);
        %fid = fopen([runpath, 'profiles/n.dat'],'w');
        %fprintf(fid, '%.15e\t%.15e\n', ne');
        %fclose(fid);

        %change temp profile
        Te = load([runpath, 'profiles/Te.dat']);
        Te(:, 2) = Te(:, 2) .* Tfac(p);
        %fid = fopen([runpath, 'profiles/Te.dat'],'w');
        %fprintf(fid, '%.15e\t%.15e\n', Te');
        %fclose(fid);
        
        disp([studyname, ' nfac = ', num2str(nfac(o)), ' Tfac = ', num2str(Tfac(p))]);

        %BALANCE CODE
        bal = Balance(runpath, shot, time, libKiLCA);
        bal.FLAG_FORCE_PROFILES = false;
        bal.setCoil(pfile);
        bal.setEqui(gfile, fluxdatapath);
        bal.setProfiles();
        bal.setModes(m, n);
        bal.setOptions();
        if(exist([runpath, 'flre/linear-data/'], 'dir') ~= 7)
            bal.write();
                bal.run();
        end
	
        bal.loadOutput();
        bal.basicPost();
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

        figure('units', 'normalized', 'outerposition', [0, 0, 1, 1]);
        axis tight
        bal.plotProfiles();
        print(['~/Balance_Results/', studyname, filename, '_Prof'], '-dpng', '-r200')
        %print(['~/Balance_Results/', studyname, filename, '_Prof'], '-dsvg')

        figure('units', 'normalized', 'outerposition', [0, 0, 1, 1]);
        axis tight
        bal.plotDe22Res();
        print(['~/Balance_Results/', studyname, filename, '_De22Res'], '-dpng', '-r200')
        print(['~/Balance_Results/', studyname, filename, '_De22Res'], '-dsvg')

        close all
    end
end
