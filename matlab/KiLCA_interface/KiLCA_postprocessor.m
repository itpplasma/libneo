classdef KiLCA_postprocessor < KiLCA_prototype_output
%classdef KiLCA_interface
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
%##########################################################################

    %author:   Philipp Ulbl
    %created:  04.11.2019
    %modified: 04.11.2019

    properties
        
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
        
        mode    %modenumber (m, n)
        
        k_th    %wavevector in theta (r-dependent)
        k_z     %wavevector in z (scalar)
        
        %------------------------------------------------------------------
        
        Br      %complex radial magnetic field
        Bth     %complex theta magnetic field
        Bz      %complex z magnetic field
        B       %absolute magnetic field
        
        dBr     %r-derivative of complex radial magnetic field
        dBth    %r-derivative of complex theta magnetic field
        dBz     %r-derivative of complex z magnetic field
        
        %------------------------------------------------------------------
        
        Jr      %complex radial current
        Jth     %complex theta current
        Jz      %complex z current
        
        J       %absolute current
        Jpar    %parallel current
        Jperp   %perpendicular current
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
            
            %calculate currents for linear data
            obj.Calc_Currents();
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
        function plotJcyl(obj, type, varargin)
            %##############################################################
            %function plotJcyl(obj, type, varargin)
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
            plot_triple(obj, 'Jr', 'Jth', 'Jz', 'J / statAmp cm^{-2}', type, varargin{:});
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
            plot_triple(obj, 'J', 'Jpar', 'Jperp', 'J / statAmp cm^{-2}', type, varargin{:});
        end
    end
    
    methods (Access = private)
        
        function calc_Currents(obj)
            %calculates current densities
            
            %save needed radius vector and parameters
            obj.r = obj.ldata.R;
            obj.rmin = min(obj.r);
            obj.rres = obj.ldata.res;
            obj.rpl = obj.parent.background.rpl;
            obj.ra = obj.parent.antenna.ra;
            obj.rmax = max(obj.r);
            obj.Rtor = obj.parent.background.Rtor;
            
            %modenumber
            obj.mode = obj.ldata.mode;
            %wavevectors
            obj.k_th = obj.mode(1) ./ obj.r;
            obj.k_z = obj.mode(2) / obj.Rtor;
            
            %extract complex magnetic field components
            obj.Br  = obj.ldata.Br_Re(:, 2)  + 1i .* obj.ldata.Br_Im(:, 2);
            obj.Bth = obj.ldata.Bth_Re(:, 2) + 1i .* obj.ldata.Bth_Im(:, 2);
            obj.Bz  = obj.ldata.Bz_Re(:, 2)  + 1i .* obj.ldata.Bz_Im(:, 2);
            %calculate absolute magnetic field
            obj.B = sqrt(obj.Br.^2 + obj.Bth.^2 + obj.Bz.^2);
            
            %calculate r-derivative of B components
            obj.dBr = gradient(obj.Br, obj.r);
            obj.dBth = gradient(obj.Bth, obj.r);
            obj.dBz = gradient(obj.Bz, obj.r);

            %calculate current densities in cylindrical coordinates
            obj.Jr  = 1i/ (4*pi) .* (obj.k_th .* obj.Bz - obj.k_z .* obj.Bth);
            obj.Jth = 1 / (4*pi) .* (1i .* obj.k_z .* obj.Br - obj.dBz);
            obj.Jz  = 1 / (4*pi) .* (obj.Bth ./ obj.r + obj.dBth - 1i .* obj.k_th .* obj.Br);

            %calculate absolute, parallel and perpendicular current
            obj.J = sqrt(obj.Jr.^2 + obj.Jth.^2 + obj.Jz.^2);
            obj.Jpar = (obj.Jr .* obj.Br + obj.Jth .* obj.Bth + obj.Jz .* obj.Bz) ./ obj.B;
            obj.Jperp = obj.J - obj.Jpar;
        end
    end
end

