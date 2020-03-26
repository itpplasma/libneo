classdef profile_preprocessor < handle & hdf5_output
%##########################################################################
%classdef profile_preprocessor < handle & hdf5_output
%##########################################################################
% description of class:
%--------------------------------------------------------------------------
% This class implements the full profile processing routines. Requires
% equilibrium to be preprocessed for preprocessing of raw profiles. Can
% also load already existing processed profiles.
%##########################################################################
% properties:
%--------------------------------------------------------------------------
% *) LIB_BALANCE, FLAG_FORCE_NEO2, FLAG_USE_NEVILLE, NEVILLE_DEG
% *) SMOOTH_LEVEL, SMOOTH_FAC, SMOOTH_METHOD
% *) Te_inf, Ti_inf, ne_inf_rel, vt_inf_rel, vth_inf, Jth_inf, Jz_inf
% READONLY:
% *) profpath, gfile, pfile, convexfile, fluxdatapath
% *) neprof, Teprof, Tiprof, vtprof, kprof
% *) b_tor, r_big, equilrqpsi, r, q, psi_pol_norm
% *) r_min, r_max, r_out
% *) Bth_out, Bz_out, B_out, pi_out, dpi_out, dTi_out, k, vhat
% *) Te, Ti, ne, vt, vth, Er, qp, Jz, Jth
%##########################################################################
% methods:
%--------------------------------------------------------------------------
% *) obj = profile_preprocessor(ppath, libbal)
% *) set_equilibrium(obj, gf, pf, cf, fpath)
% *) set_profiles(obj, ne, Te, Ti, vt)
% *) process(obj, rmin, rmax)
%########################################################################## 

%author:   Philipp Ulbl
%created:  08.01.2020

    properties
        FLAG_FORCE_NEO2 = false;            %forces recalculation of k-profile by NEO2 (within calculation of profiles)

        FLAG_USE_NEVILLE = false;           %use neville polynomials instead of splines
        NEVILLE_DEG = 9;                    %degree of neville polynomials

        SMOOTH_LEVEL = 2;                   %level of smoothing (# of derivatives smoothed)
        SMOOTH_FAC = 1e-2;                  %factor for smoothing (how much % of data is smoothed)
        SMOOTH_METHOD = 'loess';            %method of smoothing -> see MATLAB doc of smooth
        
        Te_inf = 50;        %value of Te at the edge
        Ti_inf = 50;        %value of Ti at the edge
        ne_inf_rel = 1e-1;  %value of ne at the edge relative to last
        vt_inf_rel = 1e-3;  %value of vt at the edge relative to last
        vth_inf = 0;        %value of vth at the edge
        
        Jth_inf = 0; %no current in vacuum
        Jz_inf = 0;  %no current in vacuum
    end
    properties(SetAccess = public)
        LIB_BALANCE     %path to libbalance
    
        profpath        %path to profile output
        gfile           %path to gfile
        pfile           %path to coil file
        convexfile      %path to convex file
        fluxdatapath    %path to fluxdata
        
        neprof  %path to experimental ne profile
        Teprof  %path to experimental Te profile
        Tiprof  %path to experimental Ti profile
        vtprof  %path to experimental vt profile
        kprof   %path to k profile
        
        b_tor   %toroidal magnetic field at center
        r_big   %large torus radius
        
        equilrqpsi      %data from file equil r q psi in fluxdatapath
        r               %equivalent radius
        q               %safety factor
        psi_pol_norm    %normalized poloidal flux
        
        r_min       %minimum radius for profiles
        r_max       %maximum radius for profiles
        r_out       %output radius vector
        
        Bth_out     %poloidal magnetic field used for q
        Bz_out      %toroidal magnetic field used for q
        B_out       %total magnetic field used for q
        
        pi_out       %total pressure
        dpi_out      %derivative of total pressure
        dTi_out     %derivative of ion temperature
        
        k    %k profile vector (on r)
        vhat %vhat (part of vth)
        
        Te	%profile class for Te
        Ti	%profile class for Ti          
        ne	%profile class for ne
        vt	%profile class for vt
        Jth	%profile class for Jth
        Jz	%profile class for Jz
        qp	%profile class for qp
        vth	%profile class for vth
        Er	%profile class for Er
        
    end
    properties(Constant)
        
        %##################################################################
        % CONSTANTS IN CGS
        %##################################################################

        c = 29979245800.0;      % sol in cm/s
        kB = 1.3807e-16;        % Boltzmann constant
        eVK = 1.1604e4;         % eV -> deg(K)
        echarge = 4.8031e-10;   % electron charge in statcoulomb
        emass = 9.1094e-28;     % electron mass in g
        pmass = 1.6726e-24;     % proton mass in g
        eV = 1.6022e-12;        % electron volts in erg    
    end
    
    methods
        function obj = profile_preprocessor(ppath)
            %##############################################################
            %function obj = profile_preprocessor(ppath)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % constructor of the profile_preprocessor class.
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % ppath ... path of profiles
            %############################################################## 
            
            obj.profpath = ppath;
            
            %get path to LIB BALANCE which is 1 directory above this file
            s = which(mfilename);      %filename
            s = [fileparts(s), '/../'];%directory above
            s = what(s);               %get absolute path
            obj.LIB_BALANCE = [s.path, '/'];
        end
        
        function set_equilibrium(obj, gf, pf, cf, fpath)
            %##############################################################
            %function set_equilibrium(obj, gf, pf, cf, fpath)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % sets equilibrium information needed for profile preprocessing
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % gf    ... location of gfile
            % pf    ... location of pfile (coil file from Kisslinger)
            % cf    ... location of convexfile
            % fpath ... directory containing fluxdata (amn.dat, etc.)
            %############################################################## 
            
            %check input files
            if(exist(gf, 'file') ~= 2), error(['gfile not found in: ', gf]); end
            if(exist(pf, 'file') ~= 2), error(['pfile not found in: ', pf]); end
            if(exist(cf, 'file') ~= 2), error(['convexfile not found in: ', cf]); end
            if(exist(fpath, 'dir') ~= 7), error(['fluxdatapath not found in: ', fpath]); end
            
            obj.gfile = gf;
            obj.pfile = pf;
            obj.convexfile = cf;
            obj.fluxdatapath = fpath;
            
            %load rbig, btor
            raw = load([obj.fluxdatapath, 'btor_rbig.dat']);
            obj.b_tor = raw(1);
            obj.r_big = raw(2);

            %load equil file
            obj.equilrqpsi = dlmread([obj.fluxdatapath, 'equil_r_q_psi.dat'], '', 3, 0);
            obj.r = obj.equilrqpsi(:, 1); %equivalent radius
            obj.q = obj.equilrqpsi(:, 2); %safety factor
            obj.psi_pol_norm = obj.equilrqpsi(:, 3)./obj.equilrqpsi(end, 3); %normalized poloidal flux

        end
        
        function set_profiles(obj, ne, Te, Ti, vt)
            %##############################################################
            %function set_profiles(obj, ne, Te, Ti, vt)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % sets experimental profiles
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % ne    ... location of ne profile (on rho pol from experiment)
            % Te    ... location of Te profile (on rho pol from experiment)
            % Ti    ... location of Ti profile (on rho pol from experiment)
            % vt    ... location of vt profile (on rho pol from experiment)
            %############################################################## 
            
            %check input profiles
            if(exist(ne, 'file') ~= 2), error(['neprof not found in: ', ne]); end
            if(exist(Te, 'file') ~= 2), error(['Teprof not found in: ', Te]); end
            if(exist(Ti, 'file') ~= 2), error(['Tiprof not found in: ', Ti]); end
            if(exist(vt, 'file') ~= 2), error(['vtprof not found in: ', vt]); end
            
            obj.neprof = ne;
            obj.Teprof = Te;
            obj.Tiprof = Ti;
            obj.vtprof = vt;
            
        end
        
        function process(obj, rmin, rmax)
            %##############################################################
            %function process(obj, rmin, rmax)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % preprocess profiles and store in profile path
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % rmin  ... minimum radius for profile
            % rmax  ... maximum radius for profile
            %############################################################## 
            
            obj.r_max = rmax;
            obj.r_min = rmin;
            
            disp(['Start of Profile Preprocessor at ', datestr(datetime)])
    
            %#1: preprocess experimental profiles ne, Te, Ti, vt
            obj.process_nTv();
            %#2: preprocess safety factor q
            obj.process_q();
            %#3: preprocess neo-2 k and Er, vth
            obj.process_kEv();
            
            obj.write();
            
            disp(['Finished Profiles at ', datestr(datetime)])
        end
        
        function loadExisting(obj)
            %##############################################################
            %function loadExisting(obj)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % loads existing profiles from output path.
            %##############################################################

            obj.Te = profile_magician('Te', '', obj.profpath, 1);
            obj.Te.loadExisting();
            
            obj.Ti = profile_magician('Ti', '', obj.profpath, 1);
            obj.Ti.loadExisting();
            
            obj.ne = profile_magician('n', '', obj.profpath, 1);
            obj.ne.loadExisting();
            
            obj.vt = profile_magician('Vz', '', obj.profpath, 1);
            obj.vt.loadExisting();
            
            obj.qp = profile_magician('q', '', obj.profpath, 1);
            obj.qp.loadExisting();
            
            obj.vth = profile_magician('Vth', '', obj.profpath, 1);
            obj.vth.loadExisting();
            
            obj.Er = profile_magician('Er', '', obj.profpath, 1);
            obj.Er.loadExisting();
            
            obj.r_out = obj.Te.r_out;
        end
    end
    
    methods(Access = private)
        
        function process_nTv(obj)

            %electron temperature Te
            obj.Te = profile_magician('Te', obj.Teprof, obj.profpath, 1);
            obj.Te.read(true);
            obj.Te.d = 0.25;
            obj.Te.process(obj.r, obj.psi_pol_norm, obj.r_min, obj.r_max, obj.Te_inf);
            
            %density n
            obj.ne = profile_magician('n', obj.neprof, obj.profpath, 1e-6);
            obj.ne.read(true);
            ne_inf = obj.ne.y_in(end) * obj.ne_inf_rel; %last value from IDA is fine
            obj.ne.d = 0.3;
            obj.ne.dr_cut = 1;
            obj.ne.process(obj.r, obj.psi_pol_norm, obj.r_min, obj.r_max, ne_inf);

            %ion temperature Ti
            obj.Ti = profile_magician('Ti', obj.Tiprof, obj.profpath, 1);
            obj.Ti.read(true);
            obj.Ti.d = 0.5;
            obj.Ti.dr_cut = 1;
            obj.Ti.process(obj.r, obj.psi_pol_norm, obj.r_min, obj.r_max, obj.Ti_inf);

            %toroidal velocity vt
            obj.vt = profile_magician('Vz', obj.vtprof, obj.profpath, obj.r_big);
            obj.vt.read(true);
            obj.vt.d = 0.5; %wider cut region
            vt_inf = obj.vt.y_in(end) * obj.vt_inf_rel; %last value
            obj.vt.process(obj.r, obj.psi_pol_norm, obj.r_min, obj.r_max, vt_inf);

            %radius for all output profiles
            obj.r_out = obj.Te.r_out;
            
            %smoothing of experimental data
            quant = {'Te', 'ne', 'Ti', 'vt'};
            for j = 1:numel(quant)
                smo = smooth2level(obj.(quant{j}).y_out, obj.r_out, obj.SMOOTH_LEVEL, obj.SMOOTH_FAC, obj.SMOOTH_METHOD);
                obj.(quant{j}).y_out = smo{1};
            end
                
        end
        
        function process_q(obj)
            
            %1) calc B out of q; J out of B
            %2) do profile magic on J
            %3) recalc B out of J; recalc q out of B
            r_ode = obj.r;
            psi_ode = obj.psi_pol_norm;
            q_ode = -obj.q;

            Te_ode = interp1(obj.Te.r_out, obj.Te.y_out, r_ode, 'spline');
            Ti_ode = interp1(obj.Ti.r_out, obj.Ti.y_out, r_ode, 'spline');
            ne_ode = interp1(obj.ne.r_out, obj.ne.y_out, r_ode, 'spline');

            %total pressure
            p_ode = ne_ode .* (Te_ode + Ti_ode) .* obj.kB .* obj.eVK;
            %g-factor
            g_ode = 1 + (r_ode ./ -q_ode ./ obj.r_big).^2;
            %radial derivative of pressure
            if(obj.FLAG_USE_NEVILLE == true)
                dp_ode = eval_neville_polynom_for_array(length(r_ode), r_ode, p_ode, obj.NEVILLE_DEG, r_ode, 1, 1);
            else %numerical gradient
                dp_ode = gradient(p_ode, r_ode);
            end

            %anonymous functions for ode solver
            fq  = @(x) interp1(r_ode, q_ode, x, 'spline');
            fg  = @(x) interp1(r_ode, g_ode, x, 'spline');
            fdp = @(x) interp1(r_ode, dp_ode, x, 'spline');
            %solve ode for fields (=u)
            odefun = @(x, y) -2.0 * x .* y ./ fq(x).^2 ./ fg(x) / obj.r_big^2 - 8.0 * pi * fdp(x);
            odeopt = odeset ('RelTol', 1e-12);
            [r_ode, u_ode] = ode113(odefun, r_ode, obj.b_tor^2 * g_ode(1), odeopt);
            %get field components out for u
            Bz_ode = sign(obj.b_tor) * sqrt(u_ode ./ g_ode);
            Bth_ode = r_ode .* Bz_ode ./ q_ode / obj.r_big;

            %radial derivatives of field components
            if(obj.FLAG_USE_NEVILLE == true)
                dBth_ode = eval_neville_polynom_for_array(length(r_ode), r_ode, Bth_ode, obj.NEVILLE_DEG, r_ode, 1, 1);
                dBz_ode = eval_neville_polynom_for_array(length(r_ode), r_ode, Bz_ode, obj.NEVILLE_DEG, r_ode, 1, 1);
            else %spline (num gradient pretty bad for B)
                Bth_spline = spline(r_ode, Bth_ode);
                dBth_spline = splineder(Bth_spline);
                dBth_ode = ppval(dBth_spline, r_ode);

                Bz_spline = spline(r_ode, Bz_ode);
                dBz_spline = splineder(Bz_spline);
                dBz_ode = ppval(dBz_spline, r_ode);
            end
            %currents
            Jth_ode = - obj.c / 4 / pi * dBz_ode;
            Jz_ode = obj.c / 4 / pi * (Bth_ode ./ r_ode + dBth_ode);

            %poloidal current Jth
            obj.Jth = profile_magician('Jth', '', obj.profpath, 1);
            obj.Jth.psi_in = psi_ode;
            obj.Jth.y_in = Jth_ode;
            obj.Jth.d = 0.2;
            obj.Jth.dr_cut = -0.3;
            obj.Jth.process(r_ode, psi_ode, 0, obj.r_max, obj.Jth_inf, 'ee');
            
            %toroidal current Jz
            obj.Jz = profile_magician('Jz', '', obj.profpath, 1);
            obj.Jz.psi_in = psi_ode;
            obj.Jz.y_in = Jz_ode;
            obj.Jz.d = 0.2;
            obj.Jz.dr_cut = -0.3;
            obj.Jz.process(r_ode, psi_ode, 0, obj.r_max, obj.Jz_inf, 'ee');

            %magnetic field
            obj.Bth_out = 4 * pi / obj.c * cumtrapz(obj.Jz.r_out, obj.Jz.r_out .* obj.Jz.y_out) ./ obj.Jz.r_out ...
                        + Bth_ode(1) .* r_ode(1) ./ obj.Jz.r_out; %initial conditions
            obj.Bz_out = -4 * pi / obj.c * cumtrapz(obj.Jth.r_out, obj.Jth.y_out) + Bz_ode(1);
            obj.B_out = sqrt(obj.Bz_out.^2 + obj.Bth_out.^2);
            
            %safety factor
            q_out = obj.r_out / obj.r_big .* obj.Bz_out ./ obj.Bth_out;
            obj.qp = profile_magician('q', '', obj.profpath, 1);
            obj.qp.r_out = obj.r_out;
            obj.qp.y_out = q_out;

        end
        
        function process_kEv(obj)
        
            %path to neo-2 run within balance dir
            kpath = [obj.profpath, 'kprof/'];

            %kprof = '/temp/heyn/BALANCE_2017/33120_5500/Neo2_kcoef/kprof.33120_5500';
            obj.kprof = [kpath, 'k.dat'];

            %run NEO2 if forced or kprof does not exist
            if(exist(obj.kprof, 'file') ~= 2 || obj.FLAG_FORCE_NEO2 == true)
                obj.run_neo2(kpath);
            end
            
            %KPROF POSTPROCESS -> V POLOIDAL, E RADIAL. kasilov2014 eq. 6!
            
            %ION (!!) pressure + gradient
            obj.pi_out = obj.ne.y_out .* obj.Ti.y_out .* obj.kB .* obj.eVK;
            obj.dpi_out = gradient(obj.pi_out, obj.r_out);
            
            %ion temperature gradient
            if(obj.FLAG_USE_NEVILLE == true)
                obj.dTi_out = eval_neville_polynom_for_array(length(obj.Ti.r_out), ...
                    obj.Ti.r_out, obj.Ti.y_out, obj.NEVILLE_DEG, obj.Ti.r_out, 1, 1) .* obj.kB .* obj.eVK;
            else
                obj.dTi_out = gradient(obj.Ti.y_out, obj.Ti.r_out) .* obj.kB .* obj.eVK;
            end
            %vth prefactor without k (kasilov2014, eq. 6)
            obj.vhat = obj.c .* obj.Bz_out .* obj.dTi_out ./ (obj.echarge .* obj.B_out.^2);

            %load k profile - output of NEO-2
            raw = load(obj.kprof);
            rneo = raw(:, 1);
            kneo = raw(:, 2);
            kneo = kneo(~isnan(rneo));
            rneo = rneo(~isnan(rneo));

            %extend k profile to rmax
            rneo(end+1) = obj.r_max;
            kneo(end+1) = 0.5.*kneo(end);
            obj.k = interp1(rneo, kneo, obj.r, 'pchip', 'extrap');

            %poloidal velocity out of k
            obj.vth = profile_magician('Vth', '', obj.profpath, 1);
            obj.vth.psi_in = obj.psi_pol_norm;
            obj.vth.y_in = obj.k .* obj.vhat(1:numel(obj.r));
            obj.vth.dr_cut = -0.2;
            obj.vth.d = 0.05;
            obj.vth.process(obj.r, obj.psi_pol_norm, obj.r_min, obj.r_max, obj.vth_inf);
            
            %radial electric field
            k_out = interp1(rneo, kneo, obj.r_out, 'pchip', 'extrap');
            Er_out = 1 ./ (obj.echarge .* obj.ne.y_out) .* obj.dpi_out - k_out .* obj.dTi_out / obj.echarge;
            %add toroidal rotation
            Er_out = Er_out - obj.r_out .* obj.vt.y_out .* obj.B_out ./ (obj.r_big .* obj.c .* obj.qp.y_out);
            obj.Er = profile_magician('Er', '', obj.profpath, 1);
            obj.Er.r_out = obj.r_out;
            obj.Er.y_out = Er_out;

        end
        
        function run_neo2(obj, kpath)
            
            system(['mkdir -p ', kpath]);
            [~, ~] = system(['rm -r ', kpath, '*']);

            % PREPARE PROFILES FOR NEO-2

            %set points for neo-2
            r_border = [obj.r_min, 10;
                        11, 53;
                        58, max(obj.r)-1]; %borders -> denser in core and edge
            n_points = [5, 3, 9];   %points between borders

            %read equil r q psi
            R_eff  = obj.equilrqpsi(:, 1);
            S      = obj.equilrqpsi(:, 4) ./ obj.equilrqpsi(end, 4);
            R_beg  = obj.equilrqpsi(:, 8);
            Z_beg  = obj.equilrqpsi(:, 9);

            % generate distribution in r_eff
            r_neo = zeros(sum(n_points), 1);
            k2 = 0;
            for l = 1:size(r_border, 1)
                k1 = k2 + 1;
                k2 = k1 + n_points(l) - 1;
                r_neo(k1:k2) = linspace(r_border(l, 1), r_border(l, 2), n_points(l));
            end

            %interpolate quantities on r_eff
            s_neo  = interp1(R_eff, S, r_neo, 'pchip');
            r_beg  = interp1(R_eff, R_beg, r_neo, 'pchip');
            z_beg  = interp1(R_eff, Z_beg, r_neo, 'pchip');
            Ti_neo  = interp1(obj.Ti.r_out, obj.Ti.y_out, r_neo, 'pchip');
            ne_neo  = interp1(obj.ne.r_out, obj.ne.y_out, r_neo, 'pchip');

            %calculate kappa from ions - from plasma_params.m of Gernot Kapper
            collog  = 39.1 - 1.15 * log10(1e6 * ne_neo) + 2.3 * log10(1e-3 * Ti_neo); %coulomb logarithm - Form of Sergei
            v_th    = sqrt(2 * Ti_neo * obj.eV / obj.pmass); %thermal velocity
            tau     = 3 * obj.pmass^2 * v_th.^3 ./ (16 * sqrt(pi) * ne_neo * obj.echarge^4 .* collog); %collision time
            kappa   = -2 ./ (v_th .* tau); %kappa

            %output matrix
            M = [s_neo, r_neo, r_beg, z_beg, Ti_neo, ne_neo, kappa];
            % Write file for NEO-2 inputs
            save([kpath, 'surfaces.dat'], 'M', '-ascii')

            % PREPARE RUN FOR NEO-2

            %copy content of neo2 in libbalance to balance run path
            system(['rsync -av ', obj.LIB_BALANCE, 'neo2/', ' ', kpath]);
            %copy convexfile to kpath -> on temp (condor does not work with files
            %on proj/plasma)
            system(['cp ', obj.convexfile, ' ', kpath]);
            %create divB0 input file
            fdb0 = field_divB0(obj.gfile, obj.pfile, [kpath, 'convexwall.dat'], obj.fluxdatapath);
            fdb0.ipert = 0;
            fdb0.write([obj.LIB_BALANCE, 'blueprints/'], [kpath, 'TEMPLATE_DIR/']);

            %run python script create_surf_realspace.py
            cd(kpath);
            [~, ~] = system('python create_surf_realspace.py');

            % RUN FOR NEO-2
            run_remote_neo2([obj.LIB_BALANCE, 'neo2/remote_run.conf'], kpath); %this is way better
            %run_condor_neo2(kpath);

            % EXTRACT kprof
            abort = collect_kprof(kpath, obj.kprof);
            if(abort > 0)
                warning([num2str(abort), ' jobs not finished.'])
            end
        end
        
        function write(obj)  
            
            %write profiles to file
            quant = {'Te', 'ne', 'Ti', 'vt', 'qp', 'Er', 'vth'};
            for j = 1:numel(quant)
                obj.(quant{j}).write();
            end
        end
    end
end