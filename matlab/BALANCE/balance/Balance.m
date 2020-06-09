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
% *) EXEC_NAME, CONVEX_PATH, LIB_BALANCE, LIB_KiLCA, LIB_GPEC
% *) FLAG_RUN_TIME_EVOLUTION
% *) FLAG_FORCE_KISSLINGER, FLAG_FORCE_FOURIER, FLAG_FORCE_PROFILES
% *) FLAG_FORCE_NEO2, FLAG_FORCE_TMHD
% *) r_sep, r_ant, r_idw, r_sta, r_min, r_max, rb_cut_out, re_cut_out
% *) KiLCA_I0, KiLCA_flab, KiLCA_sigma, KiLCA_vgalsys
% *) kil_flre, kil_vacuum
% READONLY:
% *) path_run, path_factors, path_output, path_fluxdata,
%    path_profiles
% *) name, shot, time, m, n
% *) file_coil_raw, file_coil, file_equi, file_ne_raw, file_Te_raw, 
%    file_Ti_raw, file_vt_raw
% *) b_tor, r_big, r_sep_real, fdb0
% *) profiles, options, factors
% *) legacy_Atheta_antenna
% *) outputs
% *) r_res, d_res, q_res, I_res, De22_res
% *) legacy_specfactors, legacy_I_rescaled, scalefactors_sq, De22_res_resc
%##########################################################################
% methods:
%--------------------------------------------------------------------------
% *) obj = Balance(runpath, shot, time, name)
% *) setCoil(obj, pfile, cfile)
% *) setEqui(obj, gfile, fluxdatapath)
% *) setProfiles(obj, neprof, Teprof, Tiprof, vtprof)
% *) setOptions(obj, opt)
% *) setFactors(obj, f) -> not implemented in run
% *) setModes(obj, m, n)
% *) setTMHD(obj, ttype, tpath)
% *) setKiLCA(obj)
% *) write(obj)
% *) run(obj)
% *) loadOutput(obj)
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

        LIB_BALANCE     %path to libbalance
        LIB_KiLCA       %path to libkilca
        LIB_GPEC        %path to libgpec
        
        FLAG_FORCE_KISSLINGER = false;      %forces recalculation of field.dat if required
        FLAG_FORCE_FOURIER = false;         %forces recalculation of amn.dat if required
        FLAG_FORCE_PROFILES = false;        %forces recalculation of profiles
        FLAG_FORCE_NEO2 = false;            %forces recalculation of neo-2 k profile when profiles are recalculated
        FLAG_FORCE_TMHD = false;            %forces recalculation of resonant currents with toroidal mhd code
        
        FLAG_REQUIRE_KISSLINGER = false;    %code will need output of Kisslinger
        FLAG_REQUIRE_AMN = false;           %code will need output of fouriermodes with mode = 1 (amn.dat)
        
        RES_WARNING_TRIGGER = 1e-1;         %trigger value to give a warning for fluid resonances near surface resonance (distance in cm)
        
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
        KiLCA_vgalsys = -1e8;   %kilca velocity of moving frame
    end
    
    properties(SetAccess = private)
    
        path_run        %path of run of the balance code
        path_factors    %path where factor files are created
        path_output     %path where output is created
        
        path_fluxdata   %path to folder containing fluxdata
        path_profiles   %path to folder containing profiles
        path_tmhd       %path to folder containing data of toroidal mhd code
        
        name            %name of the run
        shot            %shot number
        time            %time in ms
        
        m               %mode number m
        n               %mode number n
        
        file_coil_raw   %location of raw coil file from experiment
        file_coil       %location of coil file from Kisslinger (3D)
        file_equi       %location of gfile
        
        file_tmhd       %location of file containing resonant currents of toroidal iMHD code
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
        
        tmhd            %object representing toroidal mhd code.
        kil_vacuum      %object of KiLCA vacuum run
        kil_flre        %object of KiLCA flre run
        
        zero_veperp      %zeros of veperp
        zero_vExB        %zeros of vExB
        dis_zero_veperp  %absolute distance for each mode from veperp zero
        dis_zero_vExB    %absolute distance for each mode from vExB zero
        
        outputs         %Output container of balance code (cell-array)
        
        r_res           %location of resonances for all modes
        d_res           %thickness of resonant layers for all modes
        q_res           %resonant safety factor for all modes
        I_res           %resonant total parallel current for all modes (integral over d_res)
        Da_res          %estimated Da at r_res for all modes
        De22_res        %rescaled De22 at r_res for all modes
        bifurcfactors   %factors for the bifurcation to occur (De22_res / Da_res)
        
        I_tmhd          %resonant currents computed by toroidal mhd code
        scalefactors_sq    %Scale factors Cmn^2
        
        %OLD QUANTITIES COMPUTED FOR CHECKS IF NEEDED
        
        legacy_Atheta_antenna  %poloidal vector potential from antenna
        legacy_specfactors     %Spectral factors Wmn
        legacy_I_rescaled      %rescaled resonant currents
        legacy_De22_unscaled   %unscaled De22 (recalculated using legacy_specfactors)
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
            system(['mkdir -p ', obj.path_run]);
            
            obj.path_factors = [runpath, 'factors/'];
            obj.path_output = [runpath, 'out/'];
            
            %get location of lib balance
            s = which(mfilename);      %filename
            s = [fileparts(s), '/../'];%directory above
            s = what(s);               %get absolute path
            obj.LIB_BALANCE = [s.path, '/'];
            obj.LIB_KiLCA = [obj.LIB_BALANCE, 'balance/KiLCA_interface/'];
            obj.LIB_GPEC = [obj.LIB_BALANCE, 'balance/GPEC_interface/'];
            
            if(~exist(obj.LIB_KiLCA, 'dir'))
                error('libKiLCA not found. Create symbolic link at location of this m file.')
            end
            
            addpath(obj.LIB_BALANCE);
            addpath(obj.LIB_KiLCA);
            addpath(obj.LIB_GPEC);
            
            obj.shot = shot;
            obj.time = time;
            
            if(nargin < 4 || isempty(name))
                name = 'yet another unnamed run';
            end
            obj.name = name;
        end
        
        %##################################################################
        % METHODS FOR RUNNING THE CODE
        %##################################################################
        
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
        
        function setCoil(obj, cfile, pfile)
            %##############################################################
            %function setCoil(obj, pfile, pfile)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % set coil file for balance run
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % cfile   ... location of raw coil file 
            %             (needed if pfile not exists for recalc)
            % pfile   ... (optional) location of Kisslinger coil file 
            %             (forces calculation if not present)
            %##############################################################    
            
            obj.file_coil_raw = cfile;
            
            %do this only if pfile is used
            if(nargin == 3)
                %save path of pfile
                obj.file_coil = pfile;
                %check if pfile exists
                if(exist(obj.file_coil, 'file') ~= 2 || obj.FLAG_FORCE_KISSLINGER == true)
                    coilpath = [obj.LIB_BALANCE, 'coil/'];
                    Kisslinger(coilpath, obj.shot, obj.file_coil_raw, obj.file_coil);
                end
            end
        end
        
        function setEqui(obj, gfile, fluxdatapath, calc_amn)
            %##############################################################
            %function setEqui(obj, gfile, fluxdatapath, use_amn)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % set equilibrium for balance run
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % gfile        ... location of g file
            % fluxdatapath ... directory where fluxdata is located
            % calc_amn     ... default true.
            %                  flag to calculate amn from existing pfile
            %                  (set in setCoil).
            %##############################################################    
            
            if(nargin < 4 || isempty(calc_amn))
                calc_amn = true;
            end
            
            obj.file_equi = gfile;
            obj.path_fluxdata = fluxdatapath;
            
            %check if amn exists. if not run with ipert=1 (use coil file)
            if(calc_amn == true && (exist([obj.path_fluxdata, 'amn.dat'], 'file') ~= 2 || obj.FLAG_FORCE_FOURIER == true))
                run_amn = 1;
                run_fourier = true;
            %check if btorrbig and equilrqpsi exists. if not run with ipert=0 (only equi)
            elseif(~exist([obj.path_fluxdata, 'btor_rbig.dat'], 'file') || ~exist([obj.path_fluxdata, 'equil_r_q_psi.dat'], 'file'))
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
        
        function setTMHDCode(obj, ttype, tpath)
            %##############################################################
            % setTMHDCode(obj, ttype, tpath)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % set the toroidal mhd code (tmhd) used to get scaling factors
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % ttype        ... type of tmhd code (supported: 'GPEC')
            % tpath        ... path to output of tmhd code
            %##############################################################    
            
            %check if run necessary
            if(obj.FLAG_FORCE_TMHD || ~exist(tpath, 'dir'))
                run_tmhd = true;
            else
                run_tmhd = false;
            end
                
            %type is GPEC
            if(strcmp(ttype, 'GPEC'))
                
                %save path
                obj.path_tmhd = tpath;
                
                %run if necessary
                if(run_tmhd==true)
                    %run gpec with matlab interface
                    obj.tmhd = GPEC_interface(tpath, obj.shot, obj.time, obj.name);
                    obj.tmhd.setEqui(obj.file_equi);
                    obj.tmhd.setCoil(obj.file_coil_raw);
                    obj.tmhd.write();
                    obj.tmhd.run();
                end
                
                %load data
                obj.file_tmhd = [obj.path_tmhd, 'gpec_profile_output_n2.nc'];
                %read safety factor
                qraw = ncread(obj.file_tmhd, 'q_rational');
                %read I res in A
                Iraw = ncread(obj.file_tmhd, 'I_res');
                
                %get resonant that match to the ones calculated by this
                Ires = Iraw(:, 1) + 1i .* Iraw(:, 2);
                %pick only m which are calculated in this class
                ind = ismember(abs(qraw), abs(obj.m./obj.n));
                Ires = abs(Ires(ind)') / 10; %convert to cgs (c=1)
                %fill Ires with nan if GPEC computed less modes than here
                ind = ~ismember(abs(obj.m./obj.n),abs(qraw));
                Ires(ind) = nan;
                Ires = reshape(Ires, size(obj.m));
                
                %save to class
                obj.I_tmhd = Ires;
            else
                warning('type of tmhd code not supported.')
                %save nan
                obj.I_tmhd = nan(size(obj.m));
            end
        end
        
        function setProfiles(obj, neprof, Teprof, Tiprof, vtprof, copy)
            %##############################################################
            %function setProfiles(obj, neprof, Teprof, Tiprof, vtprof, copy)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % set profiles for balance run. If copy is used, the first 4
            % parameters are ignored and profiles are loaded from copy.
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % neprof  ... location of ne profile from experiment
            % Teprof  ... location of Te profile from experiment
            % Tiprof  ... location of Ti profile from experiment
            % vtprof  ... location of vt profile from experiment
            % copy    ... location of processed profiles to use (instead of
            %             raw profiles + preprocess)
            %##############################################################    
            
            obj.path_profiles = [obj.path_run, 'profiles/'];
            %create path if not existent
            system(['mkdir -p ', obj.path_profiles]);

            if(nargin < 6 || isempty(copy))
                %check all profiles for existence by sophisticated method with cell arrays
                files = {'q', 'n', 'Te', 'Ti', 'Vz', 'Vth', 'Er'};
                paths = cellfun(@(x) [obj.path_profiles, x, '.dat'], files, 'UniformOutput', false);
                check = cellfun(@(x) exist(x, 'file') ~= 2, paths, 'UniformOutput', false);
                if(any(cell2mat(check)) || obj.FLAG_FORCE_PROFILES == true)
                    recalc_prof = true;
                else
                    recalc_prof = false;
                end
            else
                %copy profiles and 1level of subdir (e.g. k-profile)
                system(['rsync -a --exclude="/*/*/" ', copy, ' ', obj.path_profiles]);
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
                obj.profiles.FLAG_FORCE_NEO2 = obj.FLAG_FORCE_NEO2;
                obj.profiles.set_equilibrium(obj.file_equi, obj.CONVEX_PATH, obj.path_fluxdata);
                obj.profiles.set_profiles(obj.file_ne_raw, obj.file_Te_raw, obj.file_Ti_raw, obj.file_vt_raw);
                obj.profiles.process(obj.r_min, obj.r_max);
            else
                %load existing
                obj.profiles = profile_preprocessor(obj.path_profiles);
                obj.profiles.set_equilibrium(obj.file_equi, obj.CONVEX_PATH, obj.path_fluxdata);
                obj.profiles.loadExisting();
            end
        end
        
        function setKiLCA(obj)
            %##############################################################
            %function setKiLCA(obj)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % sets options for KiLCA and makes a pre-run to get the
            % resonant currents.
            %##############################################################    
            
            %set basic options
            obj.setKiLCAOptions();
            
            runs = {'kil_vacuum', 'kil_flre'};
            for k = 1:numel(runs)
            
                %write directory structure
                obj.(runs{k}).write();

                %pre-run of vacuum and flre
                obj.(runs{k}).run();

                %postprocess kilca
                obj.(runs{k}).post();
                
                %settings for balance run: single mode and profiles from
                %interface
                obj.(runs{k}).antenna.nmod = 1;
                obj.(runs{k}).antenna.write(obj.(runs{k}).BLUE_PATH, obj.(runs{k}).pathofrun);
                obj.(runs{k}).background.flag_recalc = -1;
                obj.(runs{k}).background.write(obj.(runs{k}).BLUE_PATH, obj.path_run);
            end
            
            %get data from kilca postprocessing to balance object
            obj.getKiLCAPost();
            
            %used when Amn is removed from balance
            obj.scalefactors_sq = (obj.I_tmhd./obj.I_res).^2;
            
            %rescale current
            obj.rescaleI();
            %used when Amn is included in balance
            %obj.scalefactors_sq = (obj.I_tmhd./obj.legacy_I_rescaled).^2; 
        end
        
        function setDaEstimation(obj, path_da)
            %##############################################################
            %function setDaEstimation(obj, path_da)
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
            
            %extract da at r_res
            obj.Da_res = interp1(obj.da_est.r, obj.da_est.Da, obj.r_res);
            
            %write 2 file for balance
            obj.da_est.write([obj.path_profiles, 'Da.dat']);
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
                        
            if(nargin < 2 || isempty(opt))
                %BALANCE OPTIONS
                obj.options = balanceoptions(obj.kil_flre.pathofrun, obj.kil_vacuum.pathofrun);
                obj.options.Btor = obj.b_tor;
                obj.options.Rtor = obj.r_big;
                obj.options.rmin = obj.r_sta;
                obj.options.rmax = obj.r_ant;
                obj.options.rb_cut_out = obj.rb_cut_out;
                obj.options.re_cut_out = obj.re_cut_out;
            else
                obj.options = opt;
            end
            obj.options.rsep = obj.r_sep_real;
        end
        
        function setFactors(obj, f)
            %##############################################################
            %function setFactors(obj, f)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % set factors for V-shift (Not implemented yet)
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % f     ... vector of factors
            %##############################################################    
            
            obj.factors = f;
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
            [~, ~] = system(['rsync -av ', obj.LIB_BALANCE, 'template_experimental/ ', obj.path_run]);
            %[~, ~] = system(['rsync -av ', obj.LIB_BALANCE, 'template/ ', obj.path_run]);

            %FIELD DIVB0
            obj.fdb0 = field_divB0(obj.file_equi, obj.file_coil, obj.CONVEX_PATH, obj.path_fluxdata);
            obj.fdb0.write([obj.LIB_BALANCE, 'blueprints/'], obj.path_run);

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
                kilmode = KiLCA_modes(obj.m(i), obj.n(i));
                kilmode.write(obj.kil_flre.BLUE_PATH, obj.kil_flre.path);
                
                %write antenna fac
                obj.options.antenna_fac = obj.scalefactors_sq(i);
                %obj.options.antenna_fac = 1;
                obj.options.write([obj.LIB_BALANCE, 'blueprints/'], obj.path_run);
                
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
                    system(['mv -f ' outfiles{j}, ' ', mnoutpath, ' 2>/dev/null']);
                end

                %go back
                cd(currentpath)
            end
            
            obj.loadOutput();
            
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
            	obj.outputs{i} = balanceoutput(obj.m(i), obj.n(i), obj.r_res(i));
                obj.outputs{i}.loadOutput(mnoutpath);
            end
            
            %read de22
            obj.De22_res = zeros(size(obj.m));
            
            for j = 1:numel(obj.m)
                if(~obj.options.flag_run_time_evolution)
                    obj.De22_res(j) = interp1(obj.outputs{j}.fort5000.r, obj.outputs{j}.fort5000.de22, obj.r_res(j));
                else
                    obj.De22_res(j) = interp1(obj.outputs{j}.fort5000{1}.r, obj.outputs{j}.fort5000{1}.de22, obj.r_res(j));
                end
            end
            obj.legacy_De22_unscaled = obj.De22_res;
            
            %scale de22 because content in fort5000 file is not
            obj.De22_res = obj.De22_res .* obj.scalefactors_sq;
            obj.bifurcfactors = obj.De22_res ./ obj.Da_res;
            
            obj.export2HDF5(obj.path_output, obj.name);

            %NOT NEEDED ANYMORE WITH KILCA PRE-RUN
            
%             %load vacuum kilca data
%             obj.kil_vacuum.modes.set(obj.m, obj.n);
%             obj.kil_vacuum.antenna.nmod = numel(obj.m);
%             obj.kil_vacuum.runExternal();
%             obj.kil_vacuum.loadOutput();
%             
%             %load flre kilca data
%             obj.kil_flre.modes.set(obj.m, obj.n);
%             obj.kil_flre.antenna.nmod = numel(obj.m);
%             obj.kil_flre.runExternal();
%             obj.kil_flre.loadOutput();
        end
        
        %##################################################################
        % METHODS FOR PLOTS
        %##################################################################
        
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
                    if(numel(obj.outputs{j}.fort5000) > 1)
                        f5k = obj.outputs{j}.fort5000{1};
                    else
                        f5k = obj.outputs{j}.fort5000;
                    end
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
                ph = semilogy(obj.r_res(l), obj.De22_res(l), 'o', 'MarkerSize', 5, ...
                    'HandleVisibility', 'off');
                hold on
                set(ph, 'MarkerFaceColor', get(ph, 'Color'));
                pname = ['(', num2str(obj.m(l)), ',' num2str(obj.n(l)) ,')'];
                text(obj.r_res(l) + 0.1, obj.De22_res(l), pname)
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
                    pos = get(lh, 'Position');
                    posx = 0.8;
                    posy = 0.25;
                    set(lh,'Position',[posx posy pos(3) pos(4)]);
                end
                hold off
            end
            
            set(findall(gcf,'-property','FontSize'), 'FontSize', 18)

        end
        
        %##################################################################
        % METHODS FOR EXPORT
        %##################################################################
        
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
            obj.writeHDF5(fname, '/output/', 'I_tmhd', 'resonant currents from toroidal mhd code', 'statA (c=1)');
            
            obj.writeHDF5(fname, '/output/', 'scalefactors_sq', 'scale factors Cmn^2', '1');
            obj.writeHDF5(fname, '/output/', 'De22_res', 'resonant De22', 'cm^2  s^{-1}');
            obj.writeHDF5(fname, '/output/', 'Da_res', 'resonant De22', 'cm^2  s^{-1}');
            obj.writeHDF5(fname, '/output/', 'bifurcfactors', 'De22_res / Da_res', '1');
            
            obj.writeHDF5(fname, '/output/', 'legacy_specfactors', 'spectral factors Wmn', '1');
            obj.writeHDF5(fname, '/output/', 'legacy_I_rescaled', 'resonant currents rescaled with Wmn', 'statA (c=1)');
            obj.writeHDF5(fname, '/output/', 'legacy_De22_unscaled', 'resonant De22 rescaled with Cmn', 'cm^2 s^{-1}');
            
            obj.writeHDF5(fname, '/output/', 'zero_veperp', 'zeros of veperp', 'cm');
            obj.writeHDF5(fname, '/output/', 'zero_vExB', 'zeros of vExB', 'cm');
            
            obj.writeHDF5(fname, '/output/', 'dis_zero_veperp', 'absolute distance for each mode from veperp zero', 'cm');
            obj.writeHDF5(fname, '/output/', 'dis_zero_vExB', 'absolute distance for each mode from vExB zero', 'cm');
            
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
            obj.kil_vacuum.export2HDF5(fname, '/KiLCA_vac/');
            
            %export da estimation
            obj.da_est.export2HDF5(fname, '/Da_estimation/');
            
            if(obj.options.flag_run_time_evolution)
                for k=1:numel(obj.outputs)
                    obj.outputs{k}.export2HDF5(fname, ['/time_evol/mode_m', num2str(obj.m(k)), '_n', num2str(obj.n(k)), '/']);
                end
            end
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
            
            colnames = {'SpecFac', 'I_KiLCA', 'legacy_I_rescaled', 'I_GPEC', 'FormFac', 'De22', 'De22_rescaled'};
            rownames = arrayfun(@(x,y) [num2str(x), ',', num2str(y)], obj.m, obj.n, 'UniformOutput', false);
            Igpec = sqrt(obj.legacy_specfactors) ./ obj.legacy_I_rescaled;
            tableval = {obj.legacy_specfactors', obj.I_res', obj.legacy_I_rescaled', ...
                        Igpec', obj.scalefactors_sq', obj.De22_res', obj.De22_res_resc'};

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
        
        function setKiLCAOptions(obj)
            %##############################################################
            %function setKiLCAOptions(obj, flre, vac)
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
            
            %number of modes
            nmodes = numel(obj.m);
            
            %build zones
            rad   = [obj.r_sta, obj.r_sep, obj.r_ant, obj.r_idw];
            bound = {'center', 'interface', 'antenna', 'idealwall'};
            med   = {'vacuum', 'vacuum', 'vacuum'};

            %set setting for vac and flre (only differ in run_type)
            runs = {'vacuum', 'flre'};
            for k = 1:numel(runs)
            
                runname = ['kil_', runs{k}];
                
                %standard settings of kilca
                obj.(runname) = KiLCA_interface(obj.path_run, runs{k});
                obj.(runname).BLUE_PATH = [obj.LIB_KiLCA, 'blueprints/'];
                obj.(runname).PROF_PATH = obj.path_profiles;
                obj.(runname).set_background(obj.r_big, obj.r_sep);
                obj.(runname).background.Btor = obj.b_tor;
                obj.(runname).background.flag_recalc = 1;
                obj.(runname).set_antenna(obj.r_ant, nmodes);
                obj.(runname).modes.set(obj.m, obj.n);
                med{1} = runs{k};
                obj.(runname).set_zones(rad, bound, med);

                %most common settings for balance run
                obj.(runname).background.vgalsys = obj.KiLCA_vgalsys;
                obj.(runname).antenna.I0 = obj.KiLCA_I0;
                obj.(runname).antenna.flab(1) = obj.KiLCA_flab;
                obj.(runname).zones{3}.vacuum.sigma(1) = obj.KiLCA_sigma; %resistive wall

            end
        end
        
        function getKiLCAPost(obj)
            %##############################################################
            %function getKiLCAPost(obj)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % basic postprocessing after kilca run
            %##############################################################
            
            %initialize resonant quantities
            obj.q_res = obj.m ./ obj.n;
            obj.r_res = zeros(size(obj.m));
            obj.I_res = zeros(size(obj.m));
            
            %read resonant quantities
            for j = 1:numel(obj.m)
                obj.r_res(j) = obj.kil_flre.postprocessors{j}.rres;
                obj.d_res(j) = obj.kil_flre.postprocessors{j}.d;
                obj.I_res(j) = obj.kil_flre.postprocessors{j}.Ipar;
            end
            
            obj.checkResonances();
        end
        
        function checkResonances(obj)
            %##############################################################
            %function checkResonances(obj)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % checks veperp and vExB for resonances near modes
            %##############################################################
            
            %get values from kilca postprocessor
            r = obj.kil_flre.postprocessors{1}.r;
            veperp = obj.kil_flre.postprocessors{1}.veperp;
            vExB = obj.kil_flre.postprocessors{1}.vExB;
            
            %sign of velocities
            sig_veperp = sign(veperp);
            sig_vExB = sign(vExB);
            %approximate zeros
            approx_zero_veperp = find(sig_veperp(2:end) ~= sig_veperp(1:(end-1)));
            approx_zero_vExB = find(sig_vExB(2:end) ~= sig_vExB(1:(end-1)));
            
            %initialize zero vectors
            obj.zero_veperp = nan(1, max(1, numel(approx_zero_veperp)));
            obj.zero_vExB = nan(1, max(1, numel(approx_zero_vExB)));
            obj.dis_zero_veperp = nan(numel(obj.zero_veperp), numel(obj.m));
            obj.dis_zero_vExB = nan(numel(obj.zero_veperp), numel(obj.m));
            %get zeros with fzero using approximate zeros
            %+get absolute distance for each mode
            for k = 1:numel(approx_zero_veperp)
                ind = approx_zero_veperp(k);
                obj.zero_veperp(k) = fzero(@(x) interp1(r(ind:(ind+1)), veperp(ind:(ind+1)), x), [r(ind), r(ind+1)]);
                obj.dis_zero_veperp(k,:) = abs(obj.zero_veperp(k) - obj.r_res);
            end
            for k = 1:numel(approx_zero_vExB)
                ind = approx_zero_vExB(k);
                obj.zero_vExB(k) = fzero(@(x) interp1(r(ind:(ind+1)), vExB(ind:(ind+1)), x), [r(ind), r(ind+1)]);
                obj.dis_zero_vExB(k,:) = abs(obj.zero_vExB(k) - obj.r_res);
            end
            
            %check veperp resonance
            L = obj.dis_zero_veperp < obj.RES_WARNING_TRIGGER;
            if(any(L))
                warning(['fluid resonance near mode m=', num2str(obj.m(L)),' detected']);
            end
            
            %check ExB resonance
            L = obj.dis_zero_vExB < obj.RES_WARNING_TRIGGER;
            if(any(L))
                warning(['ExB resonance near mode m=', num2str(obj.m(L)),' detected']);
            end
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
            if(~exist([obj.path_fluxdata, 'amn.dat'], 'file'))
                warning('amn.dat not found, Atheta not read.')
                obj.legacy_Atheta_antenna = nan;
                return;
            end
            
            %copy executable to read unformatted fortran file
            system(['cp ', obj.LIB_BALANCE, '/fourier/READ_amns/READ_amn_Atheta.x ', obj.path_fluxdata]);
            
            mpath = pwd();
            cd(obj.path_fluxdata);
            
            obj.legacy_Atheta_antenna = cell(1, numel(obj.m));
            
            %for each mode execute file, this creates ascii with Atheta
            for l=1:numel(obj.m)
                system(['echo ', num2str(obj.m(l)), ',', num2str(obj.n(l)), '| ./READ_amn_Atheta.x >/dev/null']);
                %load file and save data
                raw = load('A_theta.dat');
                obj.legacy_Atheta_antenna{l} = [raw(:, 1), raw(:, 3) + 1i .* raw(:, 4)];
            end
            
            cd(mpath);
        end
        
        function rescaleI(obj)
            %##############################################################
            %function rescaleI(obj)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % rescales current by spectral factor. Needs amn.dat in 
            % fluxdatapath.
            %##############################################################
            
            %read poloidal vector potential of antenna
            obj.readAntenna();
            
            %if not nan
            if(~any(cellfun(@(a) any(isnan(a), 'all'), obj.legacy_Atheta_antenna)))
                %calculate reference A from vacuum field
                Brvac = zeros(size(obj.m));
                for l=1:numel(obj.m)
                    post = obj.kil_vacuum.postprocessors{l};
                    Brvac(l) = interp1(post.r, post.Br, post.rres);
                end
                A_ref = (obj.r_res .* obj.r_big ./ obj.n) .* Brvac;
                
                %calculate realistic A (interp of Atheta antenna)
                A_realistic = zeros(size(obj.m));
                for l=1:numel(obj.m)
                    A_realistic(l) = interp1(obj.legacy_Atheta_antenna{l}(:, 1), obj.legacy_Atheta_antenna{l}(:, 2), obj.r_res(l));
                end
            %else set to nan
            else
                A_ref = nan(size(obj.I_res));
                A_realistic = nan(size(obj.I_res));
            end
            
            obj.legacy_specfactors = abs(A_realistic ./ A_ref);
            obj.legacy_I_rescaled = obj.I_res .* obj.legacy_specfactors;
        end
        
        function s = mn2string(obj, i)
            %small function to get a nice string out of m and n. i = index
             s = ['m=', num2str(obj.m(i)), ', n=', num2str(obj.n(i))];
        end
        
    end
end