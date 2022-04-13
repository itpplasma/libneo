function ProfilesRecalc(profpath, gfile, pfile, convexfile, fluxdatapath, ne, Te, Ti, vz, r_min, r_max)
%##########################################################################
%function ProfilesRecalc(profpath, gfile, pfile, convexfile, fluxdatapath, ne, Te, Ti, vz, r_min, r_max)
%##########################################################################
% description:
%--------------------------------------------------------------------------
% legacy function used for hand runs of profile preprocessing. Should be
% kept up to date with profile preprocessor class.
%##########################################################################
% input:
%--------------------------------------------------------------------------
% profpath      ... location of profiles (target location)
% gfile         ... path to gfile
% pfile         ... path to field.dat
% convexfile    ... location of convexfile
% fluxdatapath  ... fluxdatapath
% ne            ... ne profile from experiment
% Te            ... Te profile from experiment
% Ti            ... Ti profile from experiment
% vz            ... vz profile from experiment
% r_min         ... minimum radius for profiles
% r_max         ... maximum radius for profiles
%##########################################################################

    %get path to lib balance
    s = which(mfilename);      %filename
    s = [fileparts(s), '/../'];%directory above
    s = what(s);               %get absolute path
    libBalance = [s.path, '/'];
    
    FLAG_FORCE_NEO2 = false;            %forces recalculation of k-profile by NEO2 (within calculation of profiles)

    %for profile magic
    FLAG_USE_NEVILLE = false;           %use neville polynomials instead of splines
    NEVILLE_DEG = 9;                    %degree of neville polynomials

    r_max_prof = r_max + 1;
    
    %######################################################################
    % CONSTANTS IN CGS
    %######################################################################

    c = 29979245800.0;      % sol in cm/s
    kB = 1.3807e-16;        % Boltzmann constant
    eVK = 1.1604e4;         % eV -> deg(K)
    echarge = 4.8031e-10;   % electron charge in statcoulomb
    emass = 9.1094e-28;     % electron mass in g
    pmass = 1.6726e-24;     % proton mass in g
    eV = 1.6022e-12;        % electron volts in erg

    %######################################################################
    % LOAD FLUXDATA
    %######################################################################

    %output of fouriermodes: rbig, btor
    raw = load([fluxdatapath, 'btor_rbig.dat']);
    b_tor = raw(1);
    r_big = raw(2);

    %output of fouriermodes: equil r q psi
    %read equil file
    equilrqpsi = dlmread([fluxdatapath, 'equil_r_q_psi.dat'], '', 3, 0);
    r = equilrqpsi(:, 1);
    q = equilrqpsi(:, 2);
    psi_pol_norm = equilrqpsi(:, 3)./equilrqpsi(end, 3);
    
    %######################################################################
    % START
    %######################################################################

    disp(['Start of Profiles at ', datestr(datetime)])
    
    %radius for all output profiles
    r_out = Te.r_out;

    %######################################################################
    % SAFETY FACTOR
    %######################################################################

    %1) calc B out of q; J out of B
    %2) do profile magic on J
    %3) recalc B out of J; recalc q out of B
    r_ode = r;
    psi_ode = psi_pol_norm;
    q_ode = -q;

    Te_ode = interp1(Te.r_out, Te.y_out, r_ode, 'spline');
    Ti_ode = interp1(Ti.r_out, Ti.y_out, r_ode, 'spline');
    ne_ode = interp1(ne.r_out, ne.y_out, r_ode, 'spline');

    %total pressure
    p_ode = ne_ode .* (Te_ode + Ti_ode) .* kB .* eVK;
    %g-factor
    g_ode = 1 + (r_ode ./ -q_ode ./ r_big).^2;
    %radial derivative of pressure
    if(FLAG_USE_NEVILLE == true)
        dp_ode = eval_neville_polynom_for_array(length(r_ode), r_ode, p_ode, NEVILLE_DEG, r_ode, 1, 1);
    else %numerical gradient
        dp_ode = gradient(p_ode, r_ode);
        dp_ode = smooth(r_ode, dp_ode, numel(r_ode) / 200);
    end

    %anonymous functions for ode solver
    fq  = @(x) interp1(r_ode, q_ode, x, 'spline');
    fg  = @(x) interp1(r_ode, g_ode, x, 'spline');
    fdp = @(x) interp1(r_ode, dp_ode, x, 'spline');
    %solve ode for fields (=u)
    odefun = @(x, y) -2.0 * x .* y ./ fq(x).^2 ./ fg(x) / r_big^2 - 8.0 * pi * fdp(x);
    odeopt = odeset ('RelTol', 1e-12);
    [r_ode, u_ode] = ode113(odefun, r_ode, b_tor^2 * g_ode(1), odeopt);
    %get field components out for u
    Bz_ode = sign(b_tor) * sqrt(u_ode ./ g_ode);
    Bth_ode = r_ode .* Bz_ode ./ q_ode / r_big;

    %radial derivatives of field components
    if(FLAG_USE_NEVILLE == true)
        dBth_ode = eval_neville_polynom_for_array(length(r_ode), r_ode, Bth_ode, NEVILLE_DEG, r_ode, 1, 1);
        dBz_ode = eval_neville_polynom_for_array(length(r_ode), r_ode, Bz_ode, NEVILLE_DEG, r_ode, 1, 1);
    else %spline (num gradient pretty bad for B)
        Bth_spline = spline(r_ode, Bth_ode);
        dBth_spline = splineder(Bth_spline);
        dBth_ode = ppval(dBth_spline, r_ode);

        Bz_spline = spline(r_ode, Bz_ode);
        dBz_spline = splineder(Bz_spline);
        dBz_ode = ppval(dBz_spline, r_ode);
    end
    %currents
    Jth_ode = - c / 4 / pi * dBz_ode;
    Jz_ode = c / 4 / pi * (Bth_ode ./ r_ode + dBth_ode);

    %poloidal current Jth
    Jth = profile_magician('Jth', '', '', 1);
    Jth.psi_in = psi_ode;
    Jth.y_in = Jth_ode;
    Jth.d = 0.2;
    Jth.dr_cut = -0.3;
    Jth_inf = 0; %no current in vacuum
    Jth.process(r, psi_ode, 0, r_max_prof, Jth_inf, 'ee');
    %toroidal current Jz
    Jz = profile_magician('Jz', '', '', 1);
    Jz.psi_in = psi_ode;
    Jz.y_in = Jz_ode;
    Jz.d = 0.2;
    Jz.dr_cut = -0.3;
    Jz_inf = 0; %no current in vacuum
    Jz.process(r, psi_ode, 0, r_max_prof, Jz_inf, 'ee');

    %magnetic field
    Bth_out = 4 * pi / c * cumtrapz(Jz.r_out, Jz.r_out .* Jz.y_out) ./ Jz.r_out ...
                + Bth_ode(1) .* r_ode(1) ./ Jz.r_out; %initial conditions
    Bz_out = -4 * pi / c * cumtrapz(Jth.r_out, Jth.y_out) + Bz_ode(1);

    %safety factor
    q_out = r_out / r_big .* Bz_out ./ Bth_out;
    qp = profile_magician('q', '', profpath, 1);
    qp.r_out = Jz.r_out;
    qp.y_out = q_out;
    qp.write();

    %##########################################################################
    % NEO-2 -> KPROF
    %##########################################################################

    %path to neo-2 run within balance dir
    kpath = [profpath, 'kprof/'];
    
    %kprof = '/temp/heyn/BALANCE_2017/33120_5500/Neo2_kcoef/kprof.33120_5500';
    kprof = [kpath, 'k.dat'];
    
    if(exist(kprof, 'file') ~= 2 || FLAG_FORCE_NEO2 == true)
        run_neo2 = true;
    else
        run_neo2 = false;
    end
    
    %run NEO2 if true
    if(run_neo2 == true)
        
        system(['mkdir -p ', kpath]);
        [~, ~] = system(['rm -r ', kpath, '*']);
        
        % PREPARE PROFILES FOR NEO-2

        %set points for neo-2
        r_border = [r_min, 10;
                    11, 53;
                    58, max(r)-1]; %borders -> denser in core and edge
        n_points = [5, 3, 9];   %points between borders

        %read equil r q psi
        R_eff  = equilrqpsi(:, 1);
        S      = equilrqpsi(:, 4) ./ equilrqpsi(end, 4);
        R_beg  = equilrqpsi(:, 8);
        Z_beg  = equilrqpsi(:, 9);

        % generate distribution in r_eff
        r_neo = zeros(sum(n_points), 1);
        k2 = 0;
        for k = 1:size(r_border, 1)
            k1 = k2 + 1;
            k2 = k1 + n_points(k) - 1;
            r_neo(k1:k2) = linspace(r_border(k, 1), r_border(k, 2), n_points(k));
        end

        %interpolate quantities on r_eff
        s_neo  = interp1(R_eff, S, r_neo, 'pchip');
        r_beg  = interp1(R_eff, R_beg, r_neo, 'pchip');
        z_beg  = interp1(R_eff, Z_beg, r_neo, 'pchip');
        Ti_neo  = interp1(Ti.r_out, Ti.y_out, r_neo, 'pchip');
        ne_neo  = interp1(ne.r_out, ne.y_out, r_neo, 'pchip');

        %calculate kappa from ions - from plasma_params.m of Gernot Kapper
        collog  = 39.1 - 1.15 * log10(1e6 * ne_neo) + 2.3 * log10(1e-3 * Ti_neo); %coulomb logarithm - Form of Sergei
        v_th    = sqrt(2 * Ti_neo * eV / pmass); %thermal velocity
        tau     = 3 * pmass^2 * v_th.^3 ./ (16 * sqrt(pi) * ne_neo * echarge^4 .* collog); %collision time
        kappa   = -2 ./ (v_th .* tau); %kappa

        %output matrix
        M = [s_neo, r_neo, r_beg, z_beg, Ti_neo, ne_neo, kappa];
        % Write file for NEO-2 inputs
        save([kpath, 'surfaces.dat'], 'M', '-ascii')

        % PREPARE RUN FOR NEO-2

        %copy content of neo2 in libbalance to balance run path
        system(['rsync -av ', libBalance, 'neo2/', ' ', kpath]);
        %copy convexfile to kpath -> on temp (condor does not work with files
        %on proj/plasma)
        system(['cp ', convexfile, ' ', kpath]);
        %create divB0 input file
        fdb0 = field_divB0(gfile, pfile, [kpath, 'convexwall.dat'], fluxdatapath);
        fdb0.ipert = 0;
        fdb0.write([libBalance, 'blueprints/'], [kpath, 'TEMPLATE_DIR/']);

        %run python script create_surf_realspace.py
        cd(kpath);
        [~, ~] = system('python create_surf_realspace.py');

        % RUN FOR NEO-2
        run_remote_neo2([libBalance, 'neo2/remote_run.conf'], kpath); %this is way better
        %run_condor_neo2(kpath);

        % EXTRACT kprof
        abort = collect_kprof(kpath, kprof);
        if(abort > 0)
            warning([num2str(abort), ' jobs not finished.'])
        end
    end
    
    %##########################################################################
    % KPROF POSTPROCESS -> V POLOIDAL, E RADIAL
    %##########################################################################

    B_out = sqrt(Bz_out.^2 + Bth_out.^2);
    %p_out = ne.y_out .* (Te.y_out + Ti.y_out) .* kB .* eVK;
    pi_out = ne.y_out .* (Ti.y_out) .* kB .* eVK; %IONS!
    dpi_out = gradient(pi_out, r_out);
    dpi_out = smooth(r_out, dpi_out, numel(r_out) / 200);
    %ion temperature gradient
    if(FLAG_USE_NEVILLE == true)
        dTi_out = eval_neville_polynom_for_array(length(Ti.r_out), Ti.r_out, Ti.y_out, NEVILLE_DEG, Ti.r_out, 1, 1) .* kB .* eVK;
    else
        dTi_out = gradient(Ti.y_out, Ti.r_out) .* kB .* eVK;
        dTi_out = smooth(Ti.r_out, dTi_out, numel(r_out) / 200);
    end
    %vth prefactor without k (kasilov2014, eq. 6)
    vdihat = c .* Bz_out .* dTi_out ./ (echarge .* B_out.^2);

    %load k profile - output of NEO-2
    raw = load(kprof);
    rneo = raw(:, 1);
    kneo = raw(:, 2);
    kneo = kneo(~isnan(rneo));
    rneo = rneo(~isnan(rneo));

    %extend k profile to rmax
    rneo(end+1) = r_max;
    kneo(end+1) = 0.5.*kneo(end);
    kneo_interp = interp1(rneo, kneo, r, 'pchip', 'extrap');

    %poloidal velocity out of k
    vth = profile_magician('Vth', '', profpath, 1);
    vth.psi_in = psi_pol_norm;
    vth.y_in = kneo_interp .* vdihat(1:numel(r));
    vth_inf = 0;
    vth.dr_cut = -0.2;
    vth.d = 0.05;
    vth.process(r, psi_pol_norm, r_min, r_max_prof, vth_inf);
    vth.write();

    %radial electric field
    k_out = interp1(rneo, kneo, r_out, 'pchip', 'extrap');
    Er_out = 1 ./ (echarge .* ne.y_out) .* dpi_out - k_out .* dTi_out / echarge;
    %add toroidal rotation
    Er_out = Er_out - r_out .* vz.y_out .* B_out ./ (r_big .* c .* qp.y_out);
    Er = profile_magician('Er', '', profpath, 1);
    Er.r_out = r_out;
    Er.y_out = Er_out;
    Er.write();
    
    disp(['Finished Profiles at ', datestr(datetime)])
end
