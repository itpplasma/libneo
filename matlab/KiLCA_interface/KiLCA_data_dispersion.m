classdef KiLCA_data_dispersion < KiLCA_prototype_output
%classdef KiLCA_data_dispersion
%##########################################################################
% description of class:
%--------------------------------------------------------------------------
% A class purely used for storage of KiLCA dispersiondata which is a result 
% of a run, when nmod in antenna.in is greater than 0 and disp in 
% output.in = 2. The data is stored in dependent properties, which means 
% not all the data is in the memory until used.
%##########################################################################
% properties:
%--------------------------------------------------------------------------
% *) path
% *) R
%##########################################################################
% methods:
%--------------------------------------------------------------------------
% *) function obj = KiLCA_data_dispersion(path)
% *) function plot(obj, q, varargin)
%##########################################################################

%author:   Philipp Ulbl
%created:  22.08.2019
%modified: 22.08.2019
        
    properties
        path        %path of modenumber subfolder containing EB.dat
    end
    
    properties (Dependent=true)
        
        R           %radius
        
        k1_Re       %real part 1
        k1_Im       %imag part 1
        k2_Re       %real part 2
        k2_Im       %imag part 2
        k3_Re       %real part 3
        k3_Im       %imag part 3
        k4_Re       %real part 4
        k4_Im       %imag part 4
    end
    
    %get accessors of the dependent properties
    methods     
        
        %------------------------------------------------------------------
        
        function q = get.R(obj)
           q = get_quantity(obj, 1);
        end
        
        function q = get.k1_Re(obj)
           q = [obj.R, get_quantity(obj, 2)];
        end
        function q = get.k1_Im(obj)
           q = [obj.R, get_quantity(obj, 3)];
        end
        function q = get.k2_Re(obj)
           q = [obj.R, get_quantity(obj, 4)];
        end
        function q = get.k2_Im(obj)
           q = [obj.R, get_quantity(obj, 5)];
        end
        function q = get.k3_Re(obj)
           q = [obj.R, get_quantity(obj, 6)];
        end
        function q = get.k3_Im(obj)
           q = [obj.R, get_quantity(obj, 7)];
        end
        function q = get.k4_Re(obj)
           q = [obj.R, get_quantity(obj, 8)];
        end
        function q = get.k4_Im(obj)
           q = [obj.R, get_quantity(obj, 9)];
        end
    end    

    methods (Access = public)
        function obj = KiLCA_data_dispersion(path)
            %##############################################################
            %function obj = KiLCA_data_dispersion(path)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % constructor of the class. needs the path of the modenumber
            % subfolder in dispersion-data.
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % path  ... path of modenumber subfolder in dispersion-data
            %##############################################################
            
            obj.path = path;
        end
        
        function plot(obj, varargin)
            %##############################################################
            %function plot(obj, varargin)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % plots magnetic field components as 3 subplots in a row
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % q         ... quantity to plot (propertyname)
            % varargin  ... plot arguments (if used, type must be non-empty!)
            %##############################################################
            
            figure('units', 'normalized', 'outerposition', [0, 0, 1, 1]);
            
            %plot real parts
            subplot(1, 2, 1)
            hold on
            for k=1:4
                p = plot_single(obj, ['k', num2str(k), '_Re'], 'k / cm^{-1}', varargin{:});
                set(p, 'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 2, 'MarkerFaceColor', get(p,'Color')); 
            end
            hold off
            
            %plot imag parts
            subplot(1, 2, 2)
            hold on
            for k=1:4
                p = plot_single(obj, ['k', num2str(k), '_Im'], 'k / cm^{-1}', varargin{:});
                set(p, 'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 2, 'MarkerFaceColor', get(p,'Color')); 
            end
            hold off
        end
    end
    
    methods (Access = private)
        function q = get_quantity(obj, col)
            %returns the quantity of the specified column in zone_0_kr.dat
            q = importdata([obj.path, 'zone_0_kr.dat']);
            q = q(:, col);
        end
    end
end