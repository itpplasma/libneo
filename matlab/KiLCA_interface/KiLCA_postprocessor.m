classdef KiLCA_postprocessor < KiLCA_prototype_output
    %KILCA_POSTPROCESSOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        parent
        bdata
        ldata
        
        r
        rmin
        rres
        rpl
        ra
        rmax
        Rtor
        
        mode
        
        k_th
        k_z
        
        Br
        Bth
        Bz
        B
        
        dBr
        dBth
        dBz
        
        Jr
        Jth
        Jz
        
        J
        Jpar
        Jperp
    end
    
    methods (Access = public)
        function obj = KiLCA_postprocessor(parent, imode)
            %KILCA_POSTPROCESSOR Construct an instance of this class
            %   Detailed explanation goes here
            
            obj.parent = parent;
            
            obj.bdata = parent.backgrounddata;
            obj.ldata = parent.lineardata{imode};
            
            obj.Calc_Currents();
        end
                
        function Calc_Currents(obj)
            
            obj.r = obj.ldata.R;
            obj.rmin = min(obj.r);
            obj.rres = obj.ldata.res;
            obj.rpl = obj.parent.background.rpl;
            obj.ra = obj.parent.antenna.ra;
            obj.rmax = max(obj.r);
            obj.Rtor = obj.parent.background.Rtor;
            
            obj.mode = obj.ldata.mode;
            
            obj.k_th = obj.mode(1) ./ obj.r;
            obj.k_z = obj.mode(2) / obj.Rtor;
            
            obj.Br  = obj.ldata.Br_Re(:, 2)  + 1i .* obj.ldata.Br_Im(:, 2);
            obj.Bth = obj.ldata.Bth_Re(:, 2) + 1i .* obj.ldata.Bth_Im(:, 2);
            obj.Bz  = obj.ldata.Bz_Re(:, 2)  + 1i .* obj.ldata.Bz_Im(:, 2);
            obj.B = sqrt(obj.Br.^2 + obj.Bth.^2 + obj.Bz.^2);
            
            obj.dBr = gradient(obj.Br, obj.r);
            obj.dBth = gradient(obj.Bth, obj.r);
            obj.dBz = gradient(obj.Bz, obj.r);

            obj.Jr  = 1i/ (4*pi) .* (obj.k_th .* obj.Bz - obj.k_z .* obj.Bth);
            obj.Jth = 1 / (4*pi) .* (1i .* obj.k_z .* obj.Br - obj.dBz);
            obj.Jz  = 1 / (4*pi) .* (obj.Bth ./ obj.r + obj.dBth - 1i .* obj.k_th .* obj.Br);

            obj.J = sqrt(obj.Jr.^2 + obj.Jth.^2 + obj.Jz.^2);
            obj.Jpar = (obj.Jr .* obj.Br + obj.Jth .* obj.Bth + obj.Jz .* obj.Bz) ./ obj.B;
            obj.Jperp = obj.J - obj.Jpar;
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
        
    end
end

