classdef profile_magician < handle & hdf5_output
%##########################################################################
%classdef profile_magician < handle & hdf5_output
%##########################################################################
% description of class:
%--------------------------------------------------------------------------
% This class incooporates "profile magic" by Martin Heyn.
% 
% Profiles are read from file (on rho or psi as x-values) and processed to
% behave well from the separatrix to the antenna. This is done by exp_cut
% to have a smooth transition from the profile value at the separatrix to a
% constant value at the end of the profile.
%##########################################################################
% properties:
%--------------------------------------------------------------------------
% *) type, factor, path_in, path_out
% *) psi_in, y_in, r_out, y_out, rho_flag
% *) dr_cut, d
% *) cut_exp
%##########################################################################
% methods:
%--------------------------------------------------------------------------
% *) obj = profile_magician(t, fin, pout, rbig)
% *) read(obj, isrho)
% *) process(obj, r, psi, r_min, r_max_prof, y_inf)
% *) write(obj)
%########################################################################## 

%author:   Philipp Ulbl
%created:  08.01.2020

    properties
        %------------------------------------------------------------------
        
        type            %type of profile
        factor          %factor for profile
        path_in         %path of input profile (on psi or rho)
        path_out        %path of output profile (on r equiv)
        
        %------------------------------------------------------------------
        
        psi_in          %raw input profile x (psi)
        y_in            %raw input profile y
        r_out           %output profile x (r equiv)
        y_out           %output profile y
        
        rho_flag        %indicates if psi_in was squared from raw data
        
        %------------------------------------------------------------------
        
        dr_cut = 0.2    %distance to maximum r before cut (from separatrix)
        d = 0.1         %width of cut
        
        %------------------------------------------------------------------
        
        cut_exp = @(x, xc, d) 1.0./(1.0 + exp((x-xc)/d)); %cut function
        cut_ee = @(x, xc, d) exp(-exp((x-xc)/d)); %cut function 2 (current)
    end
    
    methods
        function obj = profile_magician(t, fin, pout, fac)
            %##############################################################
            %function obj = profile_magician(t, fin, pout, fac)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % constructor of the profile_magician class.
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % t     ... type of profile
            % fin   ... path + name of input file
            % pout  ... path to output file
            % fac   ... factor for profile (units or rbig for Vz or etc.)
            %############################################################## 
            
            %profile type
            obj.type = t;
            
            %rbig needed as factor for vt
            obj.factor = fac;
            
            %save paths
            obj.path_in  = fin;
            obj.path_out = [pout, obj.type, '.dat'];
        end
        
        function read(obj, isrho)
            %##############################################################
            %function read(obj, isrho)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % reads raw data from file.
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % isrho ... flag if profile is given on rho instead of psi
            %           then rho is squared to get psi. Is saved in
            %           property rho_flag.
            %############################################################## 
            
            %read profile
            raw = load(obj.path_in);
            obj.psi_in = raw(:, 1);
            obj.y_in = raw(:, 2);
            
            %if profile is given on rho, square to get psi_norm
            if(isrho == true)
                obj.psi_in = obj.psi_in.^2;
            end
            
            %save flag for later information
            obj.rho_flag = isrho;
        end
        
        function process(obj, r, psi, r_min, r_max_prof, y_inf, type)
            %##############################################################
            %function process(obj, r, psi, r_min, r_max_prof, y_inf)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % processes data (=profile magic). requires read before. This
            % method is influenced by parameters r_cut and d in obj.
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % r          ... equivalent radius from equil
            % psi        ... normalized poloidal flux from equil
            % r_min      ... minimum radius of input profiles
            % r_max_prof ... maximum radius of finished profiles
            % y_inf      ... value of profile at the end (~infinity)
            % type       ... type of cut (default = exp). 
            %                'exp' or 'ee' (=double exp).
            %############################################################## 
            
            %check
            if(isempty(obj.psi_in) || isempty(obj.y_in))
                error('data must be read before processing.');
            end
            
            %read max r
            r_end = r(end);
            
            %make output grid
            dr = r(end)-r(end-1);
            obj.r_out = [r; (r_end+dr:dr:r_max_prof)']; %continue r outside separatrix
            
            %index for equil quantities up to separatrix
            ind_equil = find(r > r_min);
            %index for input files
            ind_in = find(obj.psi_in <= 1); %smaller than separatrix in normalized psi
            
            %intermediate r and y on equil quantities
            ri = obj.r_out(ind_equil);
            yi = interp1(obj.psi_in(ind_in), obj.y_in(ind_in), psi(ind_equil), 'linear', 'extrap') .* obj.factor; %to right units

            %cut parameter for cut_exp
            cut = r_end + obj.dr_cut;

            %intermediate y, extrapolated on r_out
            y_spline = interp1(ri, yi, obj.r_out, 'pchip');
            y_inf = y_inf .* obj.factor;
            
            %select cut type
            if(nargin < 7 || isempty(type) || strcmp(type, 'exp'))
                cutfun = obj.cut_exp;
            elseif (strcmp(type, 'ee'))
                cutfun = obj.cut_ee;
            end
            
            %output y by cut
            obj.y_out = y_spline .* cutfun(obj.r_out, cut, obj.d) + ...
                        y_inf .* (1 - cutfun(obj.r_out, cut, obj.d));
        end
        
        function write(obj)
            %##############################################################
            %function write(obj)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % writes processed data to file. requires process before.
            %##############################################################
            
            %check
            if(isempty(obj.r_out) || isempty(obj.y_out))
                error('data must be processed before writing.');
            end
            
            %open file
            fid = fopen(obj.path_out,'wb');
            %write profile
            fprintf(fid,'%.15e\t%.15e\n', [obj.r_out'; obj.y_out']);
            %close file
            fclose(fid);
        end
        
        function plot(obj)
            %##############################################################
            %function plot(obj)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % plots the output profile on r.
            %##############################################################
            
            figure('units', 'normalized', 'outerposition', [0, 0, 0.65, 1]);
            plot(obj.r_out, obj.y_out, 'o-k', 'LineWidth', 1, 'MarkerSize', 3, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b')
            xlabel('r / cm')
            ylabel(obj.type)
        end
        
        function loadExisting(obj)
            %##############################################################
            %function loadExisting(obj)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % loads existing profiles from output path.
            %##############################################################

            raw = load(obj.path_out);
            obj.r_out = raw(:, 1);
            obj.y_out = raw(:, 2);
        end
        
    end
    
    methods(Access = protected)
        
        function writeHDF5(obj, fname, loc, quant, desc, unit)
            %writes a single entry to hdf5 file: overwrite superclass
            
            h5create(fname, [loc, obj.type], size(obj.(quant)), 'Datatype', class(obj.(quant)));
            h5write(fname, [loc, obj.type], obj.(quant));
            h5writeatt(fname, [loc, obj.type], 'decription', desc);
            h5writeatt(fname, [loc, obj.type], 'unit', unit);
        end
    end
end