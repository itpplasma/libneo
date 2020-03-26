classdef KiLCA_postprocessor < KiLCA_prototype_output & hdf5_output
%classdef KiLCA_postprocessor < KiLCA_prototype_output & hdf5_output
%##########################################################################
% description of class:
%--------------------------------------------------------------------------
% This class is used to postprocess linear data of any KiLCA_interface.
% Thus in the constructor it must be associated to a lineardata class in
% the KiLCA_interface by using an index. Upon construction, all necessary
% quantities are extracted from the KiLCA_interface and the r-derivatives
% and current densities are calculated.
%##########################################################################
% properties:
%--------------------------------------------------------------------------
% *) parent, bdata, ldata
% *) r, rmin, rres, rpl, ra, rmax, Rtor
% *) mode, k_th, k_z
% *) Br, Bth, Bz, B
% *) dBr, dBth, dBz
% *) Jr, Jth, Jz, J, Jpar, Jperp
%##########################################################################
% methods:
%--------------------------------------------------------------------------
% *) function obj = KiLCA_postprocessor(parent, imode)
% *) function plotB(obj, type, varargin)
% *) function plotdB(obj, type, varargin)
% *) function plotJcyl(obj, type, varargin)
% *) function plotJfield(obj, type, varargin)
% *) function plotAll(obj, quant, varargin)
% *) function plotRes(obj, varargin)
% *) function export2HDF5(obj, fname, loc)
%##########################################################################

    %author:   Philipp Ulbl
    %created:  04.11.2019
    %modified: 26.03.2020
        
    properties
        
        DEBUG = false;       %flag for debug output
        OUT_FURTH = false;  %flag to write result of furth ode check in subdir
        
        %------------------------------------------------------------------
        
        parent  %parent KiLCA_interface class
        bdata   %background data of parent
        ldata   %lineardata of parent that is associated to this class
        
        %------------------------------------------------------------------
        
        r       %small radius vector
        rmin    %minimum radius (min(r))
        rres    %location of resonant surface
        rpl     %plasma radius
        ra      %loaction of antenna
        rmax    %maximum radius (max(r))
        Rtor    %big torus radius
        
        %------------------------------------------------------------------
        
        mode    %modenumber (m, n)
        
        kth     %wavevector in theta (r-dependent)
        kz      %wavevector in z (scalar)
        kpar    %wavevector parallel to B0
        
        q       %safety factor
        
        %------------------------------------------------------------------
        
        B0th    %equilibrium theta magnetic field
        B0z     %equilibrium z magnetic field
        B0      %equilibrium absolute magnetic field
        
        %------------------------------------------------------------------
        
        dB0th   %r-derivative of theta equilibrium field
        dB0z    %r-derivative of z equilibrium field
        
        %------------------------------------------------------------------
        
        J0th    %equilibrium theta current
        J0z     %equilibrium z current
        
        %------------------------------------------------------------------
        
        dp      %pressure gradient
        
        %------------------------------------------------------------------
        
        Er      %complex radial electric field
        Eth     %complex theta electric field
        Ez      %complex z electric field
        
        %------------------------------------------------------------------
        
        Br      %complex radial magnetic field
        Bth     %complex theta magnetic field
        Bz      %complex z magnetic field
        
        dBr     %r-derivative of complex radial magnetic field
        dBth    %r-derivative of complex theta magnetic field
        dBz     %r-derivative of complex z magnetic field
        
        %------------------------------------------------------------------
        
        Jr      %complex radial current
        Jth     %complex theta current
        Jz      %complex z current
        
        Jpar    %parallel current
        Jperp   %perpendicular current
        
        %------------------------------------------------------------------
        
        residual%residual of furths equation
        
        %------------------------------------------------------------------
        
        d = nan    %width of the resonant layer
        Ipar = nan %total parallel current
    end
    
    properties (Dependent=true)
        k       %magnitude of wavevector
        B       %absolute magnetic field
        E       %absolute electric field
        J       %absolute current density
        
        hth     %parallel unit vector th component
        hz      %parallel unit vector z component
        h       %absolute of parallel unit vector
    end
    
    methods
        function q = get.k(obj)
           q = sqrt(obj.kth.^2 + obj.kz.^2);
        end 
        function q = get.B(obj)
           q = sqrt(obj.Br.*conj(obj.Br) + obj.Bth.*conj(obj.Bth) + obj.Bz.*conj(obj.Bz));
        end 
        function q = get.E(obj)
            q = sqrt(obj.Er.*conj(obj.Er) + obj.Eth.*conj(obj.Eth) + obj.Ez.*conj(obj.Ez));
        end 
        function q = get.J(obj)
            q = sqrt(obj.Jr.*conj(obj.Jr) + obj.Jth.*conj(obj.Jth) + obj.Jz.*conj(obj.Jz));
        end 
        function q = get.hth(obj)
            q = obj.B0th./obj.B0;
        end 
        function q = get.hz(obj)
            q = obj.B0z./obj.B0;
        end 
        function q = get.h(obj)
            q = sqrt(obj.hth.^2 + obj.hz.^2);
        end 
    end
    
    methods (Access = public)
        function obj = KiLCA_postprocessor(parent, imode)
            %##############################################################
            %function obj = KiLCA_postprocessor(parent, imode)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % Constructor of the KiLCA_postprocessor class. Extracts all
            % needed data from parent and calculates derivatives of the
            % magnetic field components and current densities.
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % parent... KiLCA_interface class that owns this entity
            % imode ... index of modenumber to use for postprocessing
            %           (index of parent.lineardata)
            %##############################################################
            
            if(numel(imode) > 1)
                error('imode must be a scalar number.');
            end
            
            %set parent
            obj.parent = parent;
            %set background and linear data
            obj.bdata = parent.backgrounddata;
            obj.ldata = parent.lineardata{imode};
            
            %intiialize quantities in this class
            obj.initialize();
        end
        
        function plot_single(obj, a, u, type, varargin)
            %##############################################################
            %function p = plot_single(obj, a, u, type, varargin)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % plots single property over radius given by name a
            % --- adds to superclass method
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % a         ... property to be plot
            % u         ... ylabel as text
            % type      ... type of plot: Re, Im or Abs (=empty)
            % varargin  ... plot arguments (if used, type must be non-empty!)
            %##############################################################
            
            plot_single@KiLCA_prototype_output(obj, a, u, type, varargin{:});
         
            if(obj.rres ~= 0)
                %plot location of resonant layer
                hold on
                %this has to be done twice to get ylim correct: no idea why
                plot([obj.rres, obj.rres], ylim, ':m', 'LineWidth', 2, 'HandleVisibility','off');
                plot([obj.rres, obj.rres], ylim, ':m', 'LineWidth', 2, 'DisplayName', 'q = m/n');
                hold off
            end
        end
        
        function plotB(obj, type, varargin)
            %##############################################################
            %function plotB(obj, type, varargin)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % plots magnetic field components
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % type      ... type of plot: Re, Im or Abs (=empty)
            % varargin  ... plot arguments (if used, type must be non-empty!)
            %##############################################################
            
            %check type if varargin used
            if nargin > 2 && isempty(type)
                error('type must be set if varargin is used.')
            end
            
            if(nargin < 2 || isempty(type))
                type = 'Abs';
            end
            
            figure('units', 'normalized', 'outerposition', [0, 0, 1, 1]);
            plot_triple(obj, 'Br', 'Bth', 'Bz', 'B / G', type, varargin{:});
        end
        function plotdB(obj, type, varargin)
            %##############################################################
            %function plotdB(obj, type, varargin)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % plots r-derivatives of magnetic field
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % type      ... type of plot: Re, Im or Abs (=empty)
            % varargin  ... plot arguments (if used, type must be non-empty!)
            %##############################################################
            
            %check type if varargin used
            if nargin > 2 && isempty(type)
                error('type must be set if varargin is used.')
            end
            
            if(nargin < 2 || isempty(type))
                type = 'Abs';
            end
            
            figure('units', 'normalized', 'outerposition', [0, 0, 1, 1]);
            plot_triple(obj, 'dBr', 'dBth', 'dBz', 'dB/dr / statAmp cm^{-2} ?', type, varargin{:});
        end
        function plotJ(obj, type, varargin)
            %##############################################################
            %function plotJ(obj, type, varargin)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % plots current density components Jr, Jth, Jz in a single plot
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % type      ... type of plot: Re, Im or Abs (=empty)
            % varargin  ... plot arguments (if used, type must be non-empty!)
            %##############################################################
            
            %check type if varargin used
            if nargin > 2 && isempty(type)
                error('type must be set if varargin is used.')
            end
            
            if(nargin < 2 || isempty(type))
                type = 'Abs';
            end
            
            figure('units', 'normalized', 'outerposition', [0, 0, 1, 1]);
            plot_triple(obj, 'Jr', 'Jth', 'Jz', 'J / statA cm^{-2}', type, varargin{:});
        end
        function plotJfield(obj, type, varargin)
            %##############################################################
            %function plotJfield(obj, type, varargin)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % plots current density components J0, Jpar, Jperp in a single 
            % plot
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % type      ... type of plot: Re, Im or Abs (=empty)
            % varargin  ... plot arguments (if used, type must be non-empty!)
            %##############################################################
             
            %check type if varargin used
            if nargin > 2 && isempty(type)
                error('type must be set if varargin is used.')
            end
            
            if(nargin < 2 || isempty(type))
                type = 'Abs';
            end
            
            figure('units', 'normalized', 'outerposition', [0, 0, 1, 1]);
            plot_triple(obj, 'J', 'Jpar', 'Jperp', 'J / statA cm^{-2}', type, varargin{:});
        end
        
        function plotAll(obj, quant, varargin)
            %##############################################################
            %function plotAll(obj, quant, varargin)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % plots all parts of all components of given quantity in 3x3
            % plot
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % quant     ... quantity to plot: B, dB, J, E
            % varargin  ... plot arguments (if used, type must be non-empty!)
            %##############################################################
            
            %choose what components to plot and labels
            if(strcmp(quant, 'Jfield'))
                quant = 'J';
                u = 'J / statAmp cm^{-2}';
                c = {'par', 'perp', ''};
            else
                if(strcmp(quant, 'B'))
                    u = 'B / G';
                elseif(strcmp(quant, 'dB'))
                    u = 'dB/dr / statAmp cm^{-2} ?';
                elseif(strcmp(quant, 'J'))
                    u = 'J / statAmp cm^{-2}';
                elseif(strcmp(quant, 'E'))
                    u = 'E / statV cm^{-1}';
                end
                c = {'r', 'th', 'z'};
            end
            
            figure('units', 'normalized', 'outerposition', [0, 0, 1, 1]);
            %Re Br
            subplot(3, 3, 1)
            plot_single(obj, [quant, c{1}], u, 'Re', varargin{:});
            %Im Br
            subplot(3, 3, 2)
            plot_single(obj, [quant, c{1}], u, 'Im', varargin{:});
            %Abs Br
            subplot(3, 3, 3)
            plot_single(obj, [quant, c{1}], u, 'Abs', varargin{:});
            
            %Re Bth
            subplot(3, 3, 4)
            plot_single(obj, [quant, c{2}], u, 'Re', varargin{:});
            %Im Bth
            subplot(3, 3, 5)
            plot_single(obj, [quant, c{2}], u, 'Im', varargin{:});
            %Abs Bth
            subplot(3, 3, 6)
            plot_single(obj, [quant, c{2}], u, 'Abs', varargin{:});
            
            %Re Bth
            subplot(3, 3, 7)
            plot_single(obj, [quant, c{3}], u, 'Re', varargin{:});
            %Im Bth
            subplot(3, 3, 8)
            plot_single(obj, [quant, c{3}], u, 'Im', varargin{:});
            %Abs Bth
            subplot(3, 3, 9)
            plot_single(obj, [quant, c{3}], u, 'Abs', varargin{:});
            
        end
        
        function plotRes(obj, type, varargin)
            %##############################################################
            %function plotRes(obj, type, varargin)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % plots residual of furths equation
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % type      ... type of plot: Re, Im or Abs (=empty)
            % varargin  ... plot arguments (if used, type must be non-empty!)
            %##############################################################
            
            %check type if varargin used
            if nargin > 2 && isempty(type)
                error('type must be set if varargin is used.')
            end
            
            if(nargin < 2 || isempty(type))
                type = 'Abs';
            end
            
            figure('units', 'normalized', 'outerposition', [0, 0, 1, 1]);
            axis tight
            
            plot_single(obj, 'residual', 'residual', type, varargin{:});
            title('Residual of Furths equation')
            
            if(~isempty(obj.d) && ~isnan(obj.d))
                hold on
                plot([obj.rres-obj.d/2, obj.rres-obj.d/2], ylim, ':r', 'LineWidth', 2, 'DisplayName', 'r_s - d')
                plot([obj.rres+obj.d/2, obj.rres+obj.d/2], ylim, ':r', 'LineWidth', 2, 'DisplayName', 'r_s + d')
                hold off
            end
            
            legend()
        end
        
        function export2HDF5(obj, fname, loc)
            %##############################################################
            %function export2HDF5(obj, fname, loc)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % exports most important content of this class to hdf5file.
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % fname  ... name of hdf5 file with path
            % loc    ... location of this sub-hierarchy in hdf5tree
            %##############################################################
            
            obj.writeHDF5(fname, loc, 'r', 'small radius vector', 'cm');
            obj.writeHDF5(fname, loc, 'rres', 'location of resonant surface', 'cm');
            obj.writeHDF5(fname, loc, 'mode', 'modenumber (m, n)', '1');
            
            obj.writeHDF5(fname, loc, 'Er', 'complex radial electric field', 'statV cm^{-1}');
            obj.writeHDF5(fname, loc, 'Eth', 'complex poloidal electric field', 'statV cm^{-1}');
            obj.writeHDF5(fname, loc, 'Ez', 'complex toroidal electric field', 'statV cm^{-1}');
            
            obj.writeHDF5(fname, loc, 'Br', 'complex radial magnetic field', 'G');
            obj.writeHDF5(fname, loc, 'Bth', 'complex poloidal magnetic field', 'G');
            obj.writeHDF5(fname, loc, 'Bz', 'complex toroidal magnetic field', 'G');
            
            obj.writeHDF5(fname, loc, 'dBr', 'r-derivative of complex radial magnetic field', 'G cm^{-1}');
            obj.writeHDF5(fname, loc, 'dBth', 'r-derivative of complex poloidal magnetic field', 'G cm^{-1}');
            obj.writeHDF5(fname, loc, 'dBz', 'r-derivative of complex toroidal magnetic field', 'G cm^{-1}');
        
            obj.writeHDF5(fname, loc, 'Jr', 'complex radial current', 'statA cm^{-2} c=1');
            obj.writeHDF5(fname, loc, 'Jth', 'complex poloidal current', 'statA cm^{-2} c=1');
            obj.writeHDF5(fname, loc, 'Jz', 'complex toroidal current', 'statA cm^{-2} c=1');
            obj.writeHDF5(fname, loc, 'Jpar', 'complex parallel current', 'statA cm^{-2} c=1');
            
            obj.writeHDF5(fname, loc, 'residual', 'residual of furths equation', 'G cm');
            obj.writeHDF5(fname, loc, 'd', 'width of the resonant layer', 'cm');
            obj.writeHDF5(fname, loc, 'Ipar', 'total parallel current', 'statA c=1');
        end
    end
    
    methods (Access = private)
        
        function initialize(obj)
            
            %only look at radii smaller than plasma radius
            ind = obj.ldata.R < obj.parent.background.rpl;
            %ind = 1:numel(obj.ldata.R);
            
            %save needed radius vector and parameters
            [obj.r, ind, ~] = unique(obj.ldata.R(ind));
            obj.rmin = min(obj.r);
            obj.rres = obj.ldata.res;
            obj.rpl  = obj.parent.background.rpl;
            obj.ra   = obj.parent.antenna.ra;
            obj.rmax = max(obj.r);
            obj.Rtor = obj.parent.background.Rtor;
            
            %extract equilibrium magnetic field components
            obj.B0th = interp1(obj.bdata.b0th(:, 1), obj.bdata.b0th(:, 2), obj.r, 'spline');
            obj.B0z  = interp1(obj.bdata.b0z(:, 1), obj.bdata.b0z(:, 2), obj.r, 'spline');
            obj.B0   = interp1(obj.bdata.b0(:, 1), obj.bdata.b0(:, 2), obj.r, 'spline');
            
            %calculate r derivatives of equilibrium field
            obj.dB0th = gradient(obj.B0th, obj.r);
            obj.dB0z  = gradient(obj.B0z, obj.r);
            
            %calculate equilibrium current
            obj.J0th = - 1/(4*pi) .* obj.dB0z;
            obj.J0z  =   1/(4*pi) .* 1./obj.r .* gradient((obj.r .* obj.B0th), obj.r);
            
            %calculate pressure gradient from equilibrium press balance
            obj.dp = obj.J0th .* obj.B0z - obj.J0z .* obj.B0th;
            
            %modenumber
            obj.mode = obj.ldata.mode;
            %wavevectors
            obj.kth  = obj.mode(1) ./ obj.r;
            obj.kz   = obj.mode(2) / obj.Rtor;
            obj.kpar = (obj.kth .* obj.B0th + obj.kz .* obj.B0z) ./ obj.B0;
            %safety factor
            obj.q = interp1(obj.bdata.q_i(:, 1), obj.bdata.q_i(:, 2), obj.r);
            
            %extract complex electric field components
            obj.Er  = obj.ldata.Er_Re(ind, 2)  + 1i .* obj.ldata.Er_Im(ind, 2);
            obj.Eth = obj.ldata.Eth_Re(ind, 2) + 1i .* obj.ldata.Eth_Im(ind, 2);
            obj.Ez  = obj.ldata.Ez_Re(ind, 2)  + 1i .* obj.ldata.Ez_Im(ind, 2);
            
            %extract complex magnetic field components
            obj.Br  = obj.ldata.Br_Re(ind, 2)  + 1i .* obj.ldata.Br_Im(ind, 2);
            obj.Bth = obj.ldata.Bth_Re(ind, 2) + 1i .* obj.ldata.Bth_Im(ind, 2);
            obj.Bz  = obj.ldata.Bz_Re(ind, 2)  + 1i .* obj.ldata.Bz_Im(ind, 2);
            
            %calculate r-derivative of B components
            obj.dBr  = gradient(obj.Br, obj.r);
            obj.dBth = gradient(obj.Bth, obj.r);
            obj.dBz  = gradient(obj.Bz, obj.r);

            %calculate current densities in cylindrical coordinates
            obj.Jr  = 1i/ (4*pi) .* (obj.kth .* obj.Bz - obj.kz .* obj.Bth);
            obj.Jth = 1 / (4*pi) .* (1i .* obj.kz .* obj.Br - obj.dBz);
            obj.Jz  = 1 / (4*pi) .* (obj.Bth ./ obj.r + obj.dBth - 1i .* obj.kth .* obj.Br);

            %calculate parallel and perpendicular current
            obj.Jpar = (obj.Jth .* obj.B0th + obj.Jz .* obj.B0z) ./ obj.B0;
            obj.Jperp = nan; %WRONG!!!
            
            %check furths equation
            obj.check_furth();
            
            %calculate layer with
            obj.calc_layerwidth();
            
            %calculate total parallel current
            if(~isnan(obj.d))
                obj.calc_Ipar();
            end
        end
        function calc_layerwidth(obj)
            
            %set d to nan if theres no resonance or run type is vacuum 
            if(obj.rres == 0 || strcmp(obj.parent.run_type, 'vacuum'))
                obj.d = nan;
                return;
            end

            %try to fit gaussian
            ind = obj.r < obj.rpl;
            model = @(b, r) b(2) .* (1/sqrt(2*pi*b(1))) .* exp(-(r-obj.rres).^2./(2*b(1)));
            options = statset('FunValCheck', 'off', 'RobustWgtFun', 'fair');
            
            warning ('off', 'all');
            beta = nlinfit(obj.r(ind), abs(obj.Jpar(ind)), model, [1e-3, 1, 1e-3, 1], options);
            warning ('on', 'all');
            
            obj.d = real(5 * sqrt(beta(1)));
            
            %plot for diag
%             figure('units', 'normalized', 'outerposition', [0, 0, 1, 0.6]);
%             axis tight
%             plot(obj.r(ind), abs(obj.Jpar(ind)), '.', 'MarkerSize', 11, 'DisplayName', 'parallel current')
%             hold on
%             plot(obj.r(ind), model(beta, obj.r(ind)), '-', 'LineWidth', 3, 'DisplayName', 'nonlinear gaussian fit')
%             xlim([obj.rres - 5 * obj.d, obj.rres + 5 * obj.d])
%             plot(obj.rres .* [1, 1], ylim, '--m', 'LineWidth', 2, 'DisplayName', 'r_{res}')
%             plot((obj.rres-obj.d) .* [1, 1], ylim, '-.r', 'LineWidth', 2, 'DisplayName', 'r_{res} \pm d')
%             plot((obj.rres+obj.d) .* [1, 1], ylim, '-.r', 'LineWidth', 2, 'HandleVisibility', 'off')
%             xlabel('r / cm')
%             ylabel('J_{||} / statA cm^{-2} (c=1)')
%             legend()
%             title('Estimation of Resonant Layer Width')
%             set(findall(gcf,'-property','FontSize'), 'FontSize', 18)
%             print('./d_est.png', '-dpng', '-r200')
%             print('./d_est.svg', '-dsvg')
            
            %show value of d if debug is true
            if(obj.DEBUG == true)
                disp(" ")
                disp([vec2str(obj.mode, '%d'), 'rres = ', num2str(obj.rres), ...
                      ', d = ', num2str(obj.d)])
            end
        end
        function calc_Ipar(obj)
            
            %Method 1: direct integration
            ind = (obj.r <= (obj.rres + obj.d / 2)) & (obj.r >= (obj.rres - obj.d / 2));
            obj.Ipar = abs(2*pi*trapz(obj.r(ind), obj.r(ind) .* obj.Jpar(ind)));
            
            %Method 2: Sergeis estimation
            term1 = interp1(obj.r, obj.hz, obj.rres, 'spline')^2 * ...
                (interp1(obj.r, real(obj.Bth), obj.rres+obj.d/2, 'spline') - ...
                 interp1(obj.r, real(obj.Bth), obj.rres-obj.d/2, 'spline'));
            term2 = interp1(obj.r, obj.hz, obj.rres, 'spline') * ...
                    interp1(obj.r, obj.hth, obj.rres, 'spline') * ...
                   (interp1(obj.r, real(obj.Bz), obj.rres+obj.d/2, 'spline') - ...
                    interp1(obj.r, real(obj.Bz), obj.rres-obj.d/2, 'spline'));
            
            Ipar2 = abs(obj.rres .* 0.5 * (term1 - term2));
            
            %Compare both
            if(obj.DEBUG == true)
                disp([vec2str(obj.mode, '%d'), 'Ipar (integrated) = ', num2str(obj.Ipar)])
                disp([vec2str(obj.mode, '%d'), 'Ipar (estimated) = ', num2str(Ipar2)])
            end
        end
        
        function check_furth(obj)
            
            %this function checks if Furths equation is fulfilled. The
            %exact version used here is eq (5, 6) in Heyn08.
            
            %calculate Fp
            Term1 = (4*pi) .* (obj.kth .* obj.J0z - obj.kz .* obj.J0th) ./ (obj.r .* obj.k.^2);
            Term2 = obj.kz .* obj.dB0z;
            Term3 = (4*pi*obj.kz^2)./(obj.kpar.*obj.B0) .* obj.dp;
            
            Fac1 = obj.r ./ (obj.kpar.*obj.B0);
            Fac2 = 2 ./ (obj.r.^2 .* obj.k.^2);
            
            Fp = Fac1 .* (gradient(Term1, obj.r) + Fac2 .* (Term2 + Term3));
            
            %calculate left side of the equation
            EqTerm1 = gradient(obj.r .* obj.dBr ./ (obj.k.^2), obj.r);
            EqTerm2 = -obj.r.*obj.Br.*(1-gradient(1./(obj.r.*obj.k.^2)));
            EqTerm3 = -obj.r.*obj.Br.*Fp;
            %no antenna current implemented by now: will need dirac-delta
            %EqTerm4 = -4*pi.*obj.r./(1i.*obj.kth) .* obj.parent.antenna.I0;
            
            if(obj.OUT_FURTH == true)
                system('mkdir -p furth');
                data = [obj.r, real(EqTerm1), imag(EqTerm1), ...
                               real(EqTerm2), imag(EqTerm2), ...
                               real(EqTerm3), imag(EqTerm3)];
                name = ['./furth/furth_', ...
                        num2str(obj.mode(1)), '_', num2str(obj.mode(2)), ...
                        '_', obj.parent.run_type, '.dat'];
                save_file(data, name);
            end
            
            %check equation
            obj.residual = EqTerm1 + EqTerm2 + EqTerm3;% + EqTerm4;
        end
    end
end

