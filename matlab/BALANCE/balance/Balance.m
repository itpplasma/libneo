classdef Balance < handle & hdf5_output
%##########################################################################
%classdef Balance < handle & hdf5_output
%##########################################################################
% description of class:
%--------------------------------------------------------------------------
% This class is used to interact with the balance code.
%##########################################################################
% properties:
%--------------------------------------------------------------------------
% *) EXEC_NAME, CONVEX_PATH, SOFT_PATH, LIB_BALANCE, LIB_KiLCA
% *) FLAG_FORCE_KISSLINGER, FLAG_FORCE_FOURIER, FLAG_FORCE_PROFILES
% *) r_sep, r_ant, r_idw, r_sta, r_min, r_max, rb_cut_out, re_cut_out
% *) KiLCA_I0, KiLCA_flab, KiLCA_sigma
% READONLY:
% *) path_run, path_factors, path_output, path_fluxdata,
%    path_profiles
% *) name, shot, time, m, n
% *) file_coil_raw, file_coil, file_equi, file_ne_raw, file_Te_raw, 
%    file_Ti_raw, file_vt_raw
% *) b_tor, r_big, r_sep_real, fdb0
% *) profiles, options, factors
% *) kil_flre, kil_vac, Atheta_antenna
% *) outputs
% *) r_res, d_res, q_res, I_res, De22_res
% *) specfactors, I_res_resc, scalefactors, De22_res_resc
%##########################################################################
% methods:
%--------------------------------------------------------------------------
% *) obj = Balance(runpath, shot, time, libkilca, name)
% *) setCoil(obj, pfile, cfile)
% *) setEqui(obj, gfile, fluxdatapath)
% *) setProfiles(obj, neprof, Teprof, Tiprof, vtprof)
% *) setOptions(obj, opt)
% *) setFactors(obj, f) -> not implemented in run
% *) setModes(obj, m, n)
% *) write(obj)
% *) run(obj)
% *) loadOutput(obj)
% *) basicPost(obj)
% *) rescaleD(obj, fname, type)
% *) rescaleI(obj)
% *) plotDe22(obj)
% *) plotDe22Single(obj, varargin)
% *) plotDe22Res(obj)
% *) plotProfiles(obj, varargin)
% *) export2HDF5(obj, path, name)
% *) export2CurTable(obj, path, name)
%########################################################################## 

%author:   Philipp Ulbl
%created:  14.01.2020

    properties
        EXEC_NAME = 'balance.x.mpif90.openmpi_x86_64'           %name of the exec file
        CONVEX_PATH = '/proj/plasma/RMP/DATA/convexwall.dat';   %location of convexfile for AUG machine
        SOFT_PATH = '/proj/plasma/soft/';                       %location of math libraries

        LIB_BALANCE     %path to libbalance
        LIB_KiLCA       %path to libkilca
        
        FLAG_FORCE_KISSLINGER = false;      %forces recalculation of field.dat
        FLAG_FORCE_FOURIER = false;         %forces recalculation of amn.dat
        FLAG_FORCE_PROFILES = false;        %forces recalculation of profiles

        r_sep = 67;             %position of separatrix
        r_ant = 70;             %position of antenna
        r_idw = 80;             %position of ideal wall

        r_sta = 3;              %minimum radius in calculations of KiLCA
        r_min = 0.1;            %minimum radius for profiles
        r_max = 70;             %maximum radius for profiles

        rb_cut_out = 67.5;      %needed in balance.in
        re_cut_out = 68;        %needed in balance.in

        KiLCA_I0 = 1e13;        %kilca antenna current
        KiLCA_flab = 1;         %kilca frequency in lab frame
        KiLCA_sigma = 1.3e16;   %kilca conductivity of resistive wall
        
    end
    
    properties(SetAccess = private)
    
        path_run        %path of run of the balance code
        path_factors    %path where factor files are created
        path_output     %path where output is created
        
        path_fluxdata   %path to folder containing fluxdata
        path_profiles   %path to folder containing profiles
        
        name            %name of the run
        shot            %shot number
        time            %time in ms
        
        m               %mode number m
        n               %mode number n
        
        file_coil_raw   %location of raw coil file from experiment
        file_coil       %location of coil file from Kisslinger (3D)
        file_equi       %location of gfile
        
        file_Ires       %location of file containing resonant currents of iMHD
        file_Da         %location of file containing Da estimation
        
        file_ne_raw     %location of raw ne profile from experiment
        file_Te_raw     %location of raw Te profile from experiment
        file_Ti_raw     %location of raw Ti profile from experiment
        file_vt_raw     %location of raw vt profile from experiment
        
        b_tor           %toroidal magnetic field at center
        r_big           %big torus radius
        r_sep_real      %location of separatrix in equilibrium
        
        profiles        %profile preprocessor
        fdb0            %field div B0 interface
        
        da_est          %da_estimator class
        
        options         %options for balance code (represents balance.in)
        factors = 0     %factors for V-shift
        
        kil_flre        %KiLCA_interface for flre
        kil_vac         %KiLCA_interface for vacuum
        
        outputs         %Output container of balance code (cell-array)
        
        Atheta_antenna  %poloidal vector potential from antenna
        
        r_res           %location of resonances for all modes
        d_res           %thickness of resonant layers for all modes
        q_res           %resonant safety factor for all modes
        I_res           %resonant total parallel current for all modes (integral over d_res)
        De22_res        %De22 at r_res for all modes
        
        specfactors     %Spectral factors Wmn
        I_res_resc      %rescaled resonant currents
        scalefactors    %Scale factors Cmn
        De22_res_resc   %rescaled De22
    end
    
    methods
        function obj = Balance(runpath, shot, time, name)
            %##############################################################
            %function obj = Balance(runpath, shot, time, name)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % constructor of the Balance class.
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % runpath   ... path of run of the balance code
            % shot      ... shot number 
            % time      ... time in ms
            % name      ... name of the run (optional)
            %##############################################################    
            
            obj.path_run = runpath;
            obj.path_factors = [runpath, 'factors/'];
            obj.path_output = [runpath, 'out/'];
            
            %get location of lib balance
            s = which(mfilename);      %filename
            s = [fileparts(s), '/../'];%directory above
            s = what(s);               %get absolute path
            obj.LIB_BALANCE = [s.path, '/'];
            obj.LIB_KiLCA = [obj.LIB_BALANCE, 'balance/KiLCA_interface/'];
            
            if(~exist(obj.LIB_KiLCA, 'dir'))
                error('libKiLCA not found. Create symbolic link at location of this m file.')
            end
            
            addpath(genpath(obj.LIB_BALANCE));
            addpath(genpath(obj.LIB_KiLCA));
            
            obj.shot = shot;
            obj.time = time;
            
            if(nargin < 5 || isempty(name))
                name = 'default';
            end
            obj.name = name;
        end
        
        function setCoil(obj, pfile, cfile)
            %##############################################################
            %function setCoil(obj, pfile, pfile)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % set coil file for balance run
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % pfile   ... location of Kisslinger coil file
            % cfile   ... location of raw coil file 
            %             (needed if pfile not exists for recalc)
            %##############################################################    
            
            obj.file_coil = pfile;
            
            %check if pfile exists
            if(exist(obj.file_coil, 'file') ~= 2 || obj.FLAG_FORCE_KISSLINGER == true)
                obj.file_coil_raw = cfile;
                coilpath = [obj.LIB_BALANCE, 'coil/'];
                Kisslinger(coilpath, obj.shot, obj.file_coil_raw, obj.file_coil);
            end
        end
        
        function setEqui(obj, gfile, fluxdatapath)
            %##############################################################
            %function setEqui(obj, gfile, fluxdatapath)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % set equilibrium for balance run
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % gfile        ... location of g file
            % fluxdatapath ... directory where fluxdata is located
            %##############################################################    
            
            obj.file_equi = gfile;
            obj.path_fluxdata = fluxdatapath;
            
            %check if amn exists. if not run with ipert=1 (use coil file)
            if(exist([obj.path_fluxdata, 'amn.dat'], 'file') ~= 2 || obj.FLAG_FORCE_FOURIER == true)
                run_amn = 1;
                run_fourier = true;
            %check if btorrbig and equilrqpsi exists. if not run with ipert=0 (only equi)
            elseif(exist([obj.path_fluxdata, 'btor_rbig.dat'], 'file') ~= 2 || exist([obj.path_fluxdata, 'equil_r_q_psi.dat'], 'file') ~= 2)
                run_amn = 0;
                run_fourier = true;
            %else no run of fourier
            else
                run_fourier = false;
            end

            if(run_fourier==true)

                %path to exec
                fourierpath = [obj.LIB_BALANCE, 'fourier/'];

                %check input
                if(exist(obj.file_equi, 'file') ~= 2), error(['gfile not found in: ', obj.file_equi]); end
                if(exist(obj.file_coil, 'file') ~= 2), error(['coil file not found in: ', obj.file_coil]); end
                if(exist(obj.CONVEX_PATH, 'file') ~= 2), error(['convexfile not found in: ', obj.CONVEX_PATH]); end

                %init divB0 class
                f0 = field_divB0(obj.file_equi, obj.file_coil, obj.CONVEX_PATH, '');
                f0.ipert = run_amn;
                f0.write([obj.LIB_BALANCE, 'blueprints/'], fourierpath);

                %run
                Fouriermodes(fourierpath, obj.path_fluxdata, f0);
            end

            %output of fouriermodes: rbig, btor
            raw = load([fluxdatapath, 'btor_rbig.dat']);
            obj.b_tor = raw(1);
            obj.r_big = raw(2);

            raw = dlmread([fluxdatapath, 'equil_r_q_psi.dat'], '', 3, 0);
            obj.r_sep_real = raw(end, 1);
        end
        
        function setProfiles(obj, neprof, Teprof, Tiprof, vtprof)
            %##############################################################
            %function setProfiles(obj, neprof, Teprof, Tiprof, vtprof)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % set profiles for balance run
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % neprof  ... location of ne profile from experiment
            % Teprof  ... location of Te profile from experiment
            % Tiprof  ... location of Ti profile from experiment
            % vtprof  ... location of vt profile from experiment
            %##############################################################    
            
            obj.path_profiles = [obj.path_run, 'profiles/'];
            %create path if not existent
            system(['mkdir -p ', obj.path_profiles]);

            %check all profiles for existence by sophisticated method with cell arrays
            files = {'q', 'n', 'Te', 'Ti', 'Vz', 'Vth', 'Er'};
            paths = cellfun(@(x) [obj.path_profiles, x, '.dat'], files, 'UniformOutput', false);
            check = cellfun(@(x) exist(x, 'file') ~= 2, paths, 'UniformOutput', false);
            if(any(cell2mat(check)) || obj.FLAG_FORCE_PROFILES == true)
                recalc_prof = true;
            else
                recalc_prof = false;
            end

            %recalc profiles if necessary
            if(recalc_prof == true)

                obj.file_ne_raw = neprof;
                obj.file_Te_raw = Teprof;
                obj.file_Ti_raw = Tiprof;
                obj.file_vt_raw = vtprof;
                
                %Profiles(profpath, gfile, pfile, convexfile, fluxdatapath, neprof, Teprof, Tiprof, vtprof, r_min, r_max);
                obj.profiles = profile_preprocessor(obj.path_profiles);
                obj.profiles.set_equilibrium(obj.file_equi, obj.file_coil, obj.CONVEX_PATH, obj.path_fluxdata);
                obj.profiles.set_profiles(obj.file_ne_raw, obj.file_Te_raw, obj.file_Ti_raw, obj.file_vt_raw);
                obj.profiles.process(obj.r_min, obj.r_max);
            else
                %load existing
                obj.profiles = profile_preprocessor(obj.path_profiles);
                obj.profiles.set_equilibrium(obj.file_equi, obj.file_coil, obj.CONVEX_PATH, obj.path_fluxdata);
                obj.profiles.loadExisting();
            end
        end
        
        function setDaEstimation(obj, path_da)
            %##############################################################
            %function setProfiles(obj, path_da)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % estimates Da based on the content of path_da. If astra data
            % is found this is used, otherwise power data and if nothing is
            % found or path is '', the universal constant is used.
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % path_da ... location of files used for Da estimation.
            %##############################################################    
            
            %force true for now
            flag_combined = true;
            
            obj.da_est = da_estimator(obj.path_profiles, obj.path_fluxdata);
            obj.da_est.loadEstimation(obj.shot, obj.time, path_da, flag_combined);   
        end
        
        function setOptions(obj, opt)
            %##############################################################
            %function setOptions(obj)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % set options for the balance code.
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % opt   ... balanceoptions object
            %##############################################################    
                        
            obj.setKiLCA();
            
            if(nargin < 2 || isempty(opt))
                %BALANCE OPTIONS
                obj.options = balanceoptions(obj.kil_flre.pathofrun, obj.kil_vac.pathofrun);
                obj.options.Btor = obj.b_tor;
                obj.options.Rtor = obj.r_big;
                obj.options.rmin = obj.r_sta;
                obj.options.rmax = obj.r_ant;
                obj.options.rsep = obj.r_sep;
                obj.options.rb_cut_out = obj.rb_cut_out;
                obj.options.re_cut_out = obj.re_cut_out;
            else
                obj.options = opt;
            end
        end
        
        function setFactors(obj, f)
            %##############################################################
            %function setFactors(obj, f)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % set factors for V-shift
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % f     ... vector of factors
            %##############################################################    
            
            obj.factors = f;
        end
        
        function setModes(obj, m, n)
            %##############################################################
            %function setModes(obj, m, n)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % set modes to calculate with the balance code.
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % m     ... mode number m
            % n     ... mode number n
            %##############################################################    
            
            %check m, n same size
            if(size(m) ~= size(n))
                error('m and n must be of same size.')
            end
            %check m scalar or vector
            if(sum(size(m) > 1) > 1)
                error('m and n must be vectors or scalars.')
            end
            
            obj.m = m;
            obj.n = n;
        end
        
        function write(obj)
            %##############################################################
            %function write(obj)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % writes all options, properties and creates dir structure
            %##############################################################
            
            %copy content from template
            system(['mkdir -p ', obj.path_run]);
            [~, ~] = system(['rsync -av ', obj.LIB_BALANCE, 'template/ ', obj.path_run]);

            %FIELD DIVB0
            obj.fdb0 = field_divB0(obj.file_equi, obj.file_coil, obj.CONVEX_PATH, obj.path_fluxdata);
            obj.fdb0.write([obj.LIB_BALANCE, 'blueprints/'], obj.path_run);

            %write link to soft path
            system(['ln -sfT ', obj.SOFT_PATH, ' ', obj.path_run, 'soft']);

            system(['mkdir -p ', obj.path_factors]);
            %write factors file for each mode (the same by now for each mode)
            for i = 1:numel(obj.m)
                facfile = [obj.path_factors, 'factors_', num2str(obj.m(i)), '_', num2str(obj.n(i)), '.in'];
                fid = fopen(facfile, 'wb');
                fprintf(fid, '%.3f\n', obj.factors');
                fclose(fid);
            end
            
            system(['mkdir -p ', obj.path_output]);
            %delete content in output path
            system(['rm -r ', obj.path_output, '* 2>/dev/null']); %supress warnings
            
            %write balance.in file (balanceoptions)
            if(isempty(obj.options))
                obj.setOptions();
            end
            obj.options.write([obj.LIB_BALANCE, 'blueprints/'], obj.path_run);
            
            %write KiLCA directory structures
            obj.kil_flre.write();
            obj.kil_vac.write();
        end
        
        function run(obj)
            %##############################################################
            %function run(obj)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % runs balance code for each mode and moves files to subdir in
            % output directory
            %##############################################################
            
            %open log file
            start_time = datetime;
            logfile = [obj.path_output, 'balance_', strrep(datestr(start_time), ' ', '_'), '.log'];
            fid = fopen(logfile, 'w');

            disp(['Start of Balance at ', datestr(start_time)])
            disp(['Shot: ', num2str(obj.shot), ', Time: ', num2str(obj.time), 'ms, Name: ', obj.name])
            for i = 1:numel(obj.m)
                
                %factors for this run
                facfile = [obj.path_factors, 'factors_', num2str(obj.m(i)), '_', num2str(obj.n(i)), '.in'];
                system(['ln -sf ', facfile, ' ', obj.path_run, 'factors.in']);
                
                %dir for output
                mnoutpath = [obj.path_output, 'f_', num2str(obj.m(i)), '_', num2str(obj.n(i)), '/'];
                system(['mkdir -p ', mnoutpath]);
                
                %change mode in modes.in file
                obj.kil_flre.modes = KiLCA_modes(obj.m(i), obj.n(i));
                obj.kil_flre.modes.write(obj.kil_flre.BLUE_PATH, obj.kil_flre.path);
                %obj.kil_flre.write();
                
                %save current path
                currentpath = pwd();
                
                %change to run path and run
                cd(obj.path_run)
                disp(['Run for mode m = ', num2str(obj.m(i)), ', n = ', num2str(obj.n(i))])
                fix = 'LD_LIBRARY_PATH=/proj/plasma/soft/math_libs/64bit/sundials-2.6.2/lib/';
                [stat, res] = system([fix, ' ./', obj.EXEC_NAME]);
                %write to log file
                fprintf(fid, '%s\n', res);
                if(stat ~= 0)
                    error(['Error in Balance. Result = ', res, '. See log file in output directory.'])
                end

                %move files to output dir
                outfiles = balanceoutput.OUTFILES;
                for j = 1:numel(outfiles)
                    system(['mv -f ' outfiles{j}, ' ', mnoutpath]);
                end

                %go back
                cd(currentpath)
            end
            
            obj.loadOutput();
            obj.basicPost();
            
            disp(['Finished Balance at ', datestr(datetime)])
            disp(['Total runtime was ', string(datetime-start_time)])
            
            %close log file
            fclose(fid);
        end
        
        function loadOutput(obj)
            %##############################################################
            %function loadOutput(obj)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % loads output files into output which is a cell array
            % of balanceout objects. 1 object for each mode.
            %##############################################################
            
            %foreach mode
            for i = 1:numel(obj.m)
                %name of out folder for this mode
                mnoutpath = [obj.path_output, 'f_', num2str(obj.m(i)), '_', num2str(obj.n(i)), '/'];
                %create balanceoutput object in cell array and load data
            	obj.outputs{i} = balanceoutput(obj.m(i), obj.n(i));
                obj.outputs{i}.loadOutput(mnoutpath);
            end
            
            %load vacuum kilca data
            obj.kil_vac.modes.set(obj.m, obj.n);
            obj.kil_vac.antenna.nmod = numel(obj.m);
            obj.kil_vac.runExternal();
            obj.kil_vac.loadOutput();
            
            %load flre kilca data
            obj.kil_flre.modes.set(obj.m, obj.n);
            obj.kil_flre.antenna.nmod = numel(obj.m);
            obj.kil_flre.runExternal();
            obj.kil_flre.loadOutput();
        end
        
        function basicPost(obj)
            %##############################################################
            %function basicPost(obj)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % basic postprocessing after balance run
            %##############################################################
            
            %postprocess kilca
            obj.kil_vac.post();
            obj.kil_flre.post();
            
            %initialize resonant quantities
            obj.q_res = obj.m ./ obj.n;
            obj.r_res = zeros(size(obj.m));
            obj.I_res = zeros(size(obj.m));
            obj.De22_res = zeros(size(obj.m));
            
            %read resonant quantities
            for j = 1:numel(obj.m)
                obj.r_res(j) = obj.kil_flre.postprocessors{j}.rres;
                obj.d_res(j) = obj.kil_flre.postprocessors{j}.d;
                obj.I_res(j) = obj.kil_flre.postprocessors{j}.Ipar;
                obj.De22_res(j) = interp1(obj.outputs{j}.fort5000.r, obj.outputs{j}.fort5000.de22, obj.r_res(j));
            end
        end
        
        function rescaleD(obj, fname, type)
            %##############################################################
            %function rescaleD(obj, fac)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % rescales de22 by scalefactor computed from ratio of iMHD
            % resonant currents and KiLCA currents.
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % fname  ... filename where Ires are located.
            % type   ... type of Ires. Supported: GPEC (nc file)
            %##############################################################
            
            %check existence
            if(exist(fname, 'file') ~= 2)
                warning('file with resonant currents not found.')
                obj.scalefactors = nan(size(obj.De22));
                obj.De22_res_resc = nan(size(obj.De22));
                return;
            end
            
            obj.file_Ires = fname;
            
            %read for type=GPEC
            if(strcmp(type, 'GPEC'))
                qraw = ncread(fname, 'q_rational'); %read safety factor
                Iraw = ncread(fname, 'I_res'); %read I res in A
                Ires = Iraw(:, 1) + 1i .* Iraw(:, 2);
                %pick only m which are calculated in this class
                ind = ismember(abs(qraw), abs(obj.m./obj.n));
                Ires = abs(Ires(ind)') / 10; %cgs
            end
            %add read for output of Patricks code
            
            %rescale by real currents (resc) or artificial
            if(any(isnan(obj.I_res_resc)))
                obj.scalefactors = (Ires./obj.I_res).^2;
            else
                obj.scalefactors = (Ires./obj.I_res_resc).^2;
            end
            obj.De22_res_resc = obj.De22_res .* obj.scalefactors;
        end
        
        function rescaleI(obj)
            %##############################################################
            %function rescaleI(obj)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % rescales current by spectral factor. Needs fluxdatapath.
            %##############################################################
            
            %read poloidal vector potential of antenna
            obj.readAntenna();
            
            %if not nan
            if(~any(cellfun(@(a) any(isnan(a), 'all'), obj.Atheta_antenna)))
                %calculate reference A from vacuum field
                Brvac = zeros(size(obj.m));
                for l=1:numel(obj.m)
                    post = obj.kil_vac.postprocessors{l};
                    Brvac(l) = interp1(post.r, post.Br, post.rres);
                end
                A_ref = (obj.r_res .* obj.r_big ./ obj.n) .* Brvac;
                
                %calculate realistic A (interp of Atheta antenna)
                A_realistic = zeros(size(obj.m));
                for l=1:numel(obj.m)
                    A_realistic(l) = interp1(obj.Atheta_antenna{l}(:, 1), obj.Atheta_antenna{l}(:, 2), obj.r_res(l));
                end
            %else set to nan
            else
                A_ref = nan(size(obj.I_res));
                A_realistic = nan(size(obj.I_res));
            end
            
            obj.specfactors = abs(A_realistic ./ A_ref);
            obj.I_res_resc = obj.I_res .* obj.specfactors;
        end
        
        function plotDe22(obj)
            %##############################################################
            %function plotDe22(obj)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % plot De22 with a predefined setup
            %##############################################################
            
            %amount of subplots
            subs = 3;
            
            %calculate indices of m for subplots
            ind_min = 1:(subs-1):numel(obj.m);
            ind_max = circshift(ind_min, -1) - 1;
            ind_max(end) = numel(obj.m);
            
            %load output of Balance
            for l = 1:subs
                subplot(subs, 1, l)
                hold on
                for j = 1:numel(obj.m)
                    f5k = obj.outputs{j}.fort5000;
                    %plot de22
                    ph = plot(f5k.r, f5k.de22, '-', 'LineWidth', 3, 'DisplayName', obj.mn2string(j));
                    %plot line at rres
                    plot([obj.r_res(j), obj.r_res(j)], ylim, '--', 'Color', get(ph, 'Color'), 'LineWidth', 2, 'HandleVisibility', 'off');
                    if(j < ind_min(l) || j > ind_max(l))
                        set(ph, 'HandleVisibility', 'off');
                    end
                end
                %plot Da
                plot(xlim, 1e4 .* [1, 1], '--k', 'LineWidth', 2, 'DisplayName', 'D_a');
                %plot separatrix on last sub
                if(l==subs)
                    plot(obj.r_sep_real .* [1, 1], ylim, '--k', 'LineWidth', 2, 'DisplayName', 'r_{sep}')
                end
                title(['shot ', num2str(obj.shot), ' @', num2str(obj.time), 'ms'])
                legend('Location', 'eastoutside')
                xlabel('r / cm')
                ylabel('D_{e22} / cm^2 s^{-1}')
                xlim([obj.r_res(ind_min(l)) - 1, obj.r_res(ind_max(l)) + 1])
                ylim([0, 1.25*max(obj.De22_res(ind_min(l):ind_max(l)))])
                hold off
            end
            set(findall(gcf,'-property','FontSize'), 'FontSize', 18)
            
        end
        
        function plotDe22Single(obj, varargin)
            %##############################################################
            %function plotDe22Single(obj)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % plot De22 with 1 subplot for each mode
            %##############################################################
            
            %amount of subplots
            subs = numel(obj.m);
            
            %calculate rows and cols
            col = ceil(sqrt(3*subs/4));
            row = ceil(subs / col);
            
            %load output of Balance
            for l = 1:subs
                subplot(row, col, l)
                hold on
                
                f5k = obj.outputs{l}.fort5000;
                %plot de22
                ph = plot(f5k.r, f5k.de22, '-', 'LineWidth', 3, 'DisplayName', obj.mn2string(l), varargin{:});

                xlim([obj.r_res(l) - obj.d_res(l), obj.r_res(l) + obj.d_res(l)])
                ylim([0, max([1.25*obj.De22_res(l), ylim])])
                
                %plot separatrix on last sub
                if(l==subs)
                    plot(obj.r_sep_real .* [1, 1], ylim, '--k', 'LineWidth', 2, 'DisplayName', 'r_{sep}', 'HandleVisibility', 'off')
                    xlim([obj.r_res(l) - obj.d_res(l), max(obj.r_res(l) + obj.d_res(l), obj.r_sep_real)])
                end
                
                %plot line at rres
                plot([obj.r_res(l), obj.r_res(l)], ylim, '--', 'Color', get(ph, 'Color'), 'LineWidth', 2, 'HandleVisibility', 'off');
                
                %plot Da
                plot(xlim, 1e4 .* [1, 1], '--k', 'LineWidth', 2, 'DisplayName', 'D_a', 'HandleVisibility', 'off');
                
                title(obj.mn2string(l))
                xlabel('r / cm')
                ylabel('D_{e22} / cm^2 s^{-1}')
                hold off
            end
            
            set(findall(gcf,'-property','FontSize'), 'FontSize', 18)
        end
        
        function plotDe22Res(obj)
            %##############################################################
            %function plotDe22Res(obj)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % plot resonant De22 on positions of resonances
            %##############################################################
            
            hold on
            yyaxis right
            files = {'n', 'Te', 'Ti'};
            styles = {'-','--','-.'};
            for l=1:numel(files)
                prof = load([obj.path_profiles, files{l}, '.dat']);
                fac = floor(log10(max(prof(:, 2))));
                plot(prof(:,1), prof(:,2) ./ (10.^fac) , styles{l}, 'LineWidth', 2, ...
                    'DisplayName', [files{l}, ' . 10^{-', num2str(fac), '}']);
            end
            xlabel('r / cm')
            ylabel('n / cm^{-3} or T / eV')
            hold off
            yyaxis left
            for l=1:numel(obj.m)
                ph = semilogy(obj.r_res(l), obj.De22_res_resc(l), 'o', 'MarkerSize', 5, ...
                    'HandleVisibility', 'off');
                hold on
                set(ph, 'MarkerFaceColor', get(ph, 'Color'));
                name = ['(', num2str(obj.m(l)), ',' num2str(obj.n(l)) ,')'];
                text(obj.r_res(l) + 0.1, obj.De22_res_resc(l), name)
            end
            set(ph, 'DisplayName', 'D_{e22}', 'HandleVisibility', 'on')
            legend()
            xlabel('r / cm')
            ylabel('D / cm^2 s^{-1}')
                        
            plot(xlim, [1e4, 1e4], '--', 'LineWidth', 2, 'DisplayName', 'D_a');
            
            legend()
            title(['shot ', num2str(obj.shot), ' @', num2str(obj.time), 'ms'])
            hold off
            set(findall(gcf,'-property','FontSize'), 'FontSize', 18)
        end
            
        function plotProfiles(obj, varargin)
            %##############################################################
            %function plotProfiles(obj)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % plot profiles
            %##############################################################
            
            files = {'q', 'n', 'Te', 'Ti', 'Vz', 'Vth', 'Er'};
            
            row = 2;
            col = ceil(numel(files)/row);
            
            for j = 1:numel(files)
                subplot(row, col, j)
                hold on
                prof = load([obj.path_profiles, files{j}, '.dat']);
                plot(prof(:,1), prof(:,2), '.r', 'LineWidth', 2, 'Handlevisibility', 'off', varargin{:});
                if(~isempty(obj.r_res))
                    for l = 1:numel(obj.m)
                        plot([obj.r_res(l), obj.r_res(l)], ylim, '--', 'DisplayName', ['res ', obj.mn2string(l)]);
                    end
                end
                title(files{j})
                xlabel('r / cm')
                ylabel(files{j})
                %put legend on 8th subplot which is free
                if(j==numel(files))
                    lh=legend();
                    pos = get(lh,'Position');
                    posx = 0.8;
                    posy = 0.25;
                    set(lh,'Position',[posx posy pos(3) pos(4)]);
                end
                hold off
            end
            
            set(findall(gcf,'-property','FontSize'), 'FontSize', 18)

        end
        
        function export2HDF5(obj, path, name)
            %##############################################################
            %function export2HDF5(obj, path, name)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % exports most relevant data to hdf5 file.
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % path  ... path of hdf5 file.
            % name  ... name of hdf5 file.
            %##############################################################
            
            fname = [path, name, '_', num2str(obj.shot), '_', num2str(obj.time), '.hdf5'];
            
            if(exist(fname, 'file'))
                delete(fname);
            end
            
            %write output quantities
            obj.writeHDF5(fname, '/output/', 'q_res', 'resonant safety factor', '1');
            obj.writeHDF5(fname, '/output/', 'r_res', 'resonant surface location', 'cm');
            obj.writeHDF5(fname, '/output/', 'd_res', 'resonant layer thickness', 'cm');
            obj.writeHDF5(fname, '/output/', 'I_res', 'resonant currents', 'statA (c=1)');
            
            obj.writeHDF5(fname, '/output/', 'De22_res', 'resonant De22', 'cm^2  s^{-1}');
            obj.writeHDF5(fname, '/output/', 'specfactors', 'spectral factors Wmn', '1');
            obj.writeHDF5(fname, '/output/', 'I_res_resc', 'resonant currents rescaled with Wmn', 'statA (c=1)');
            obj.writeHDF5(fname, '/output/', 'scalefactors', 'form factos Cmn', '1');
            obj.writeHDF5(fname, '/output/', 'De22_res_resc', 'resonant De22 rescaled with Cmn', 'cm^2 s^{-1}');
            
            %write input quantities
            obj.writeHDF5(fname, '/input/', 'name', obj.name, 'string');
            
            obj.writeHDF5(fname, '/input/', 'shot', 'shot number', '1');
            obj.writeHDF5(fname, '/input/', 'time', 'shot time', 'ms');
            obj.writeHDF5(fname, '/input/', 'm', 'poloidal mode number', '1');
            obj.writeHDF5(fname, '/input/', 'n', 'toroidal mode number', '1');
            
            obj.writeHDF5(fname, '/input/', 'b_tor', 'toroidal magnetic field at center', 'G');
            obj.writeHDF5(fname, '/input/', 'r_big', 'big torus radius', 'cm');
            obj.writeHDF5(fname, '/input/', 'r_sep_real', 'location of separatrix in equilibrium', 'cm');
            
            %export profile preprocessor data
            obj.profiles.export2HDF5(fname, '/input/');
            
            %write profiles
            obj.profiles.writeHDF5(fname, '/profiles/', 'r_out', 'output radius vector', 'cm');
            obj.profiles.Te.writeHDF5(fname, '/profiles/', 'y_out', 'electron temperature profile', 'eV');
            obj.profiles.Ti.writeHDF5(fname, '/profiles/', 'y_out', 'ion temperature profile', 'eV');
            obj.profiles.ne.writeHDF5(fname, '/profiles/', 'y_out', 'density profile', 'cm^{-3}');
            obj.profiles.vt.writeHDF5(fname, '/profiles/', 'y_out', 'toroidal velocity profile', 'cm s^{-1}');
            obj.profiles.vth.writeHDF5(fname,'/profiles/', 'y_out', 'poloidal velocity profile', 'cm s^{-1}');
            obj.profiles.qp.writeHDF5(fname, '/profiles/', 'y_out', 'safety factor profile', '1');
            obj.profiles.Er.writeHDF5(fname, '/profiles/', 'y_out', 'radial electric field profile', 'statV cm^{-1}');
        
            %export KiLCA data
            obj.kil_flre.export2HDF5(fname, '/KiLCA_flre/');
            obj.kil_vac.export2HDF5(fname, '/KiLCA_vac/');
            
            %export da estimation
            obj.da_est.export2HDF5(fname, '/Da_estimation/');
        end
        
        function export2CurTable(obj, path, name)
            %##############################################################
            %function export2CurTable(obj, path, name)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % exports resonant data to table -> "CurTable"
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % path  ... path of file.
            % name  ... name of file.
            %##############################################################
            
            colnames = {'SpecFac', 'I_KiLCA', 'I_rescaled', 'I_GPEC', 'FormFac', 'De22', 'De22_rescaled'};
            rownames = arrayfun(@(x,y) [num2str(x), ',', num2str(y)], obj.m, obj.n, 'UniformOutput', false);
            Igpec = sqrt(obj.specfactors) ./ obj.I_res_resc;
            tableval = {obj.specfactors', obj.I_res', obj.I_res_resc', ...
                        Igpec', obj.scalefactors', obj.De22_res', obj.De22_res_resc'};

            fname = [path, name, '_', num2str(obj.shot), '_', num2str(obj.time), '_CurTable.txt'];
            fid = fopen(fname, 'w');

            %construct header and print
            row = 'm,n\t';
            for j = 1:numel(colnames)
                row = [row, pad(colnames{j}, 12), '\t'];
            end
            fprintf(fid, [row, '\n']);

            %construct rows and print
            for i = 1:numel(rownames)
                row = [rownames{i}, '\t'];
                for j = 1:numel(colnames)
                    row = [row, sprintf('%e', tableval{j}(i)), '\t'];
                end
                fprintf(fid, [row, '\n']);
            end
            fclose(fid);
        end
    end
    
    methods(Access = private)
                
        function setKiLCA(obj)
            %##############################################################
            %function setKiLCA(obj, flre, vac)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % link KiLCA objects to balancecode class.
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % flre  ... KiLCA_interface class with flre run type
            % vac   ... KiLCA_interface class with vacuum run type
            %##############################################################    
            
            %FLRE, VACUUM USING KiLCA CLASS
            nmodes = 1;
            rad   = [obj.r_sta, obj.r_sep, obj.r_ant, obj.r_idw];
            bound = {'center', 'interface', 'antenna', 'idealwall'};
            med   = {'vacuum', 'vacuum', 'vacuum'};

            obj.kil_vac = KiLCA_interface(obj.path_run, 'vacuum');
            obj.kil_vac.BLUE_PATH = [obj.LIB_KiLCA, 'blueprints/'];
            obj.kil_vac.PROF_PATH = obj.path_profiles;
            obj.kil_vac.set_background(obj.r_big, obj.r_sep);
            obj.kil_vac.background.Btor = obj.b_tor;
            obj.kil_vac.background.flag_recalc = 1;%should be set to -1 in future but problem because profiles are not written by KiLCA in backgrounddata with this option
            obj.kil_vac.set_antenna(obj.r_ant, nmodes);
            med{1} = obj.kil_vac.run_type;
            obj.kil_vac.set_zones(rad, bound, med);

            obj.kil_vac.antenna.I0 = obj.KiLCA_I0;
            obj.kil_vac.antenna.flab(1) = obj.KiLCA_flab;
            obj.kil_vac.zones{3}.vacuum.sigma(1) = obj.KiLCA_sigma; %resistive wall

            obj.kil_flre = KiLCA_interface(obj.path_run, 'flre');
            obj.kil_flre.BLUE_PATH = obj.kil_vac.BLUE_PATH;
            obj.kil_flre.PROF_PATH = obj.kil_vac.PROF_PATH;
            obj.kil_flre.set_background(obj.r_big, obj.r_sep);
            obj.kil_flre.background.Btor = obj.b_tor;
            obj.kil_flre.background.flag_recalc = 1;%should be set to -1 in future but problem because profiles are not written by KiLCA in backgrounddata with this option
            obj.kil_flre.set_antenna(obj.r_ant, nmodes);
            med{1} = obj.kil_flre.run_type;
            obj.kil_flre.set_zones(rad, bound, med);

            obj.kil_flre.antenna.I0 = obj.KiLCA_I0;
            obj.kil_flre.antenna.flab(1) = obj.KiLCA_flab;
            obj.kil_flre.zones{3}.vacuum.sigma(1) = obj.KiLCA_sigma; %resistive wall

        end
        
        function readAntenna(obj)
            %##############################################################
            %function readAntenna(obj)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % reads antenna vector potential for rescaling of resonant
            % current.
            %##############################################################
            
            %check existence
            if(exist(obj.path_fluxdata, 'dir') ~= 7)
                warning('fluxdata not found, Atheta not read.')
                obj.Atheta_antenna = nan;
                return;
            end
            
            %copy executable to read unformatted fortran file
            system(['cp ', obj.LIB_BALANCE, '/fourier/READ_amns/READ_amn_Atheta.x ', obj.path_fluxdata]);
            
            mpath = pwd();
            cd(obj.path_fluxdata);
            
            obj.Atheta_antenna = cell(1, numel(obj.m));
            
            %for each mode execute file, this creates ascii with Atheta
            for l=1:numel(obj.m)
                system(['echo ', num2str(obj.m(l)), ',', num2str(obj.n(l)), '| ./READ_amn_Atheta.x >/dev/null']);
                %load file and save data
                raw = load('A_theta.dat');
                obj.Atheta_antenna{l} = [raw(:, 1), raw(:, 3) + 1i .* raw(:, 4)];
            end
            
            cd(mpath);
        end
        
        function s = mn2string(obj, i)
            %small function to get a nice string out of m and n. i = index
             s = ['m=', num2str(obj.m(i)), ', n=', num2str(obj.n(i))];
        end
        
    end
end