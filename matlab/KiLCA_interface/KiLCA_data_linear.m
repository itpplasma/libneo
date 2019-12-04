classdef KiLCA_data_linear < KiLCA_prototype_output
%classdef KiLCA_data_linear
%##########################################################################
% description of class:
%--------------------------------------------------------------------------
% A class purely used for storage of KiLCA lineardata which is a result of
% a run, when nmod in antenna.in is greater than 0. The data is stored in
% dependent properties, which means not all the data is in the memory until
% used.
%##########################################################################
% properties:
%--------------------------------------------------------------------------
% *) path
% *) mode, res, flab, fmov
% *) R
% *) Er_Re,  Er_Im,  Er_Abs
% *) Eth_Re, Eth_Im, Eth_Abs
% *) Ez_Re,  Ez_Im,  Ez_Abs
% *) Br_Re,  Br_Im,  Br_Abs
% *) Bth_Re, Bth_Im, Bth_Abs
% *) Bz_Re,  Bz_Im,  Bz_Abs
% *) Fk0, FK1, FK2, S, Ftot
%##########################################################################
% methods:
%--------------------------------------------------------------------------
% *) function obj = KiLCA_lineardata(path)
% *) function plotB(obj, type, varargin)
% *) function plotE(obj, type, varargin)
%##########################################################################

%author:   Philipp Ulbl
%created:  19.02.2019
%modified: 04.11.2019
        
    properties
        path        %path of modenumber subfolder containing EB.dat
        mode        %modenumbers
        res         %pos of res layer
        flab        %frequency of lab frame
        fmov        %frequency of moving frame
    end
    
    properties (Dependent=true)
        
        R           %radius
        
        Er_Re       %imag of E radial -> 2th col of EB.dat
        Er_Im       %real of E radial -> 3th col of EB.dat
        Er_Abs      %abs of E radial-> calculated using Er_Re, Er_Im
        
        Eth_Re      %imag of E theta -> 4th col of EB.dat
        Eth_Im      %real of E theta -> 5th col of EB.dat
        Eth_Abs     %abs of E theta-> calculated using Eth_Re, Eth_Im
        
        Ez_Re       %imag of E z -> 6th col of EB.dat
        Ez_Im       %real of E z -> 7th col of EB.dat
        Ez_Abs      %abs of E z-> calculated using Ez_Re, Ez_Im
        
        Br_Re       %imag of B radial -> 8th col of EB.dat
        Br_Im       %real of B radial -> 9th col of EB.dat
        Br_Abs      %abs of B radial-> calculated using Br_Re, Br_Im
        
        Bth_Re      %imag of B theta -> 10th col of EB.dat
        Bth_Im      %real of B theta -> 11th col of EB.dat
        Bth_Abs     %abs of B theta-> calculated using Bth_Re, Bth_Im
        
        Bz_Re       %imag of B z -> 12th col of EB.dat
        Bz_Im       %real of B z -> 13th col of EB.dat
        Bz_Abs      %abs of B z-> calculated using Bz_Re, Bz_Im
        
        Fk0
        Fk1
        Fk2
        S
        Ftot
    end
    
    %get accessors of the dependent properties
    methods
        function q = get.mode(obj)
           raw = read_in([obj.path, 'mode_data.dat']);
           raw = strsplit(raw{5});
           q = str2double(raw);
        end
        function q = get.res(obj)
           raw = read_in([obj.path, 'mode_data.dat']);
           raw = strsplit(raw{8});
           q = str2double(raw{1});
        end
        function q = get.flab(obj)
           raw = read_in([obj.path, 'mode_data.dat']);
           raw = strsplit(raw{6});
           q = str2double(raw);
        end
        function q = get.fmov(obj)
           raw = read_in([obj.path, 'mode_data.dat']);
           raw = strsplit(raw{7});
           q = str2double(raw);
        end
        
        %------------------------------------------------------------------
        
        function q = get.R(obj)
           q = get_EB_quantity(obj, 1);
        end
        
        %E radial
        function q = get.Er_Re(obj)
           q = [obj.R, get_EB_quantity(obj, 2)];
        end
        function q = get.Er_Im(obj)
           q = [obj.R, get_EB_quantity(obj, 3)];
        end
        function q = get.Er_Abs(obj)
           q = [obj.R, sqrt(obj.Er_Re(:, 2).^2 + obj.Er_Im(:, 2).^2)];
        end
        
        %------------------------------------------------------------------
        
        %E theta
        function q = get.Eth_Re(obj)
           q = [obj.R, get_EB_quantity(obj, 4)];
        end
        function q = get.Eth_Im(obj)
           q = [obj.R, get_EB_quantity(obj, 5)];
        end
        function q = get.Eth_Abs(obj)
           q = [obj.R, sqrt(obj.Eth_Re(:, 2).^2 + obj.Eth_Im(:, 2).^2)];
        end
        
        %------------------------------------------------------------------
        
        %E z
        function q = get.Ez_Re(obj)
           q = [obj.R, get_EB_quantity(obj, 6)];
        end
        function q = get.Ez_Im(obj)
           q = [obj.R, get_EB_quantity(obj, 7)];
        end
        function q = get.Ez_Abs(obj)
           q = [obj.R, sqrt(obj.Ez_Re(:, 2).^2 + obj.Ez_Im(:, 2).^2)];
        end
        
        %------------------------------------------------------------------
        
        %B radial
        function q = get.Br_Re(obj)
           q = [obj.R, get_EB_quantity(obj, 8)];
        end
        function q = get.Br_Im(obj)
           q = [obj.R, get_EB_quantity(obj, 9)];
        end
        function q = get.Br_Abs(obj)
           q = [obj.R, sqrt(obj.Br_Re(:, 2).^2 + obj.Br_Im(:, 2).^2)];
        end
        
        %------------------------------------------------------------------
        
        %B theta
        function q = get.Bth_Re(obj)
           q = [obj.R, get_EB_quantity(obj, 10)];
        end
        function q = get.Bth_Im(obj)
           q = [obj.R, get_EB_quantity(obj, 11)];
        end
        function q = get.Bth_Abs(obj)
           q = [obj.R, sqrt(obj.Bth_Re(:, 2).^2 + obj.Bth_Im(:, 2).^2)];
        end
        %------------------------------------------------------------------
        
        %B z
        function q = get.Bz_Re(obj)
           q = [obj.R, get_EB_quantity(obj, 12)];
        end
        function q = get.Bz_Im(obj)
           q = [obj.R, get_EB_quantity(obj, 13)];
        end
        function q = get.Bz_Abs(obj)
           q = [obj.R, sqrt(obj.Bz_Re(:, 2).^2 + obj.Bz_Im(:, 2).^2)];
        end
        
        %------------------------------------------------------------------
        
        %fluxes
        
        function q = get.Fk0(obj)
           q1 = importdata([obj.path, 'zone_0_kin_flux_0_e.dat']);
           q2 = importdata([obj.path, 'zone_0_kin_flux_0_i.dat']);
           q3 = importdata([obj.path, 'zone_0_kin_flux_0_t.dat']);
           q = {q1(:,1), q1(:,2), q2(:,2), q3(:,2)};
        end
        function q = get.Fk1(obj)
           q1 = importdata([obj.path, 'zone_0_kin_flux_1_e.dat']);
           q2 = importdata([obj.path, 'zone_0_kin_flux_1_i.dat']);
           q3 = importdata([obj.path, 'zone_0_kin_flux_1_t.dat']);
           q = {q1(:,1), q1(:,2), q2(:,2), q3(:,2)};
        end
        function q = get.Fk2(obj)
           q1 = importdata([obj.path, 'zone_0_kin_flux_2_e.dat']);
           q2 = importdata([obj.path, 'zone_0_kin_flux_2_i.dat']);
           q3 = importdata([obj.path, 'zone_0_kin_flux_2_t.dat']);
           q = {q1(:,1), q1(:,2), q2(:,2), q3(:,2)};
        end
        
        function q = get.S(obj)
           q = importdata([obj.path, 'zone_0_poy_flux.dat']);
           q = [q(:,1), q(:,2)];
        end
        
        function q = get.Ftot(obj)
           q = importdata([obj.path, 'zone_0_tot_flux.dat']);
           q = [q(:,1), q(:,2)];
        end
        
        %------------------------------------------------------------------
    end    

    methods (Access = public)
        function obj = KiLCA_data_linear(path)
            %##############################################################
            %function obj = KiLCA_data_linear(path)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % constructor of the class. needs the path of the modenumber
            % subfolder in linear-data that contains EB.dat.
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % path  ... path of modenumber subfolder in linear-data
            %##############################################################
                        
            obj.path = path;
        end
        
        function plot_single(obj, a, u, type, varargin)
            %adds to superclass method
            
            plot_single@KiLCA_prototype_output(obj, a, u, type, varargin{:});
         
            if(obj.res ~= 0)
                %plot location of resonant layer
                hold on
                plot([obj.res, obj.res], ylim, ':m', 'LineWidth', 2, 'DisplayName', 'q = m/n');
                hold off
            end
        end
        
        function plotB(obj, type, varargin)
            %##############################################################
            %function plotB(obj, type, varargin)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % plots magnetic field components as 3 subplots in a row
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % type      ... type of plot: Re, Im or Abs (=empty)
            % varargin  ... plot arguments (if used, type must be non-empty!)
            %##############################################################
            
            %check type if varargin used
            if nargin > 1 && isempty(type)
                error('type must be set if varargin is used.')
            end
            
            figure('units', 'normalized', 'outerposition', [0, 0, 1, 1]);
            
            if nargin > 1 && strcmp(type, 'Re')
                plot_triple(obj, 'Br_Re', 'Bth_Re', 'Bz_Re', 'B / G', type, varargin{:});
                
            elseif nargin > 1 && strcmp(type, 'Im')
                plot_triple(obj, 'Br_Im', 'Bth_Im', 'Bz_Im', 'B / G', type, varargin{:});
                
            else
                plot_triple(obj, 'Br_Abs', 'Bth_Abs', 'Bz_Abs', 'B / G', 'Abs', varargin{:});
            end
        end
        
        function plotE(obj, type, varargin)
            %##############################################################
            %function plotE(obj, type, varargin)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % plots electric field components as 3 subplots in a row
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % type      ... type of plot: Re, Im or Abs (=empty)
            % varargin  ... plot arguments (if used, type must be non-empty!)
            %##############################################################
            
            %check type if varargin used
            if nargin > 1 && isempty(type)
                error('type must be set if varargin is used.')
            end
            
            figure('units', 'normalized', 'outerposition', [0, 0, 1, 1]);
            
            if nargin > 1 && strcmp(type, 'Re')
                plot_triple(obj, 'Er_Re', 'Eth_Re', 'Ez_Re', 'E / statV cm^{-1}', type, varargin{:});
                
            elseif nargin > 1 && strcmp(type, 'Im')
                plot_triple(obj, 'Er_Im', 'Eth_Im', 'Ez_Im', 'E / statV cm^{-1}', type, varargin{:});
                
            else
                plot_triple(obj, 'Er_Abs', 'Eth_Abs', 'Ez_Abs', 'E / statV cm^{-1}', 'Abs', varargin{:});
            end
        end
    end
    
    methods (Access = private)
        function q = get_EB_quantity(obj, col)
            %returns the quantity of the specified column in EB.dat
            q = importdata([obj.path, 'EB.dat']);
            q = q(:, col);
        end
    end
end