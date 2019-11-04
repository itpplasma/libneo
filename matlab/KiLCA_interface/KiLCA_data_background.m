classdef KiLCA_data_background < KiLCA_prototype_output
%classdef KiLCA_data_background
%##########################################################################
% description of class:
%--------------------------------------------------------------------------
% A class purely used for storage of KiLCA backgroundata which is a result
% of every KiLCA run. The data is stored in dependent properties, which 
% means not all the data is in the memory until used.
%##########################################################################
% properties:
%--------------------------------------------------------------------------
% *) path
% *) b0, b0th, b0z
% *) dPHI0, Er_i, hth, hz
% *) j0p_m, j0s_m, j0th, j0th_m
% *) ne_m, ne_p, n_i, ni_m, ni_p
% *) nue, nui, q_i
% *) Te_i, Te_m, Ti_i, Ti_m
% *) Vpe_m, Vpe_p, Vpi_m, Vpi_p, Vs0e_m, Vs0i_m, Vse_m, Vsi_m, Vte_p,
% Vth_i, Vth_m, Vti_p, Vz_i, Vz_m
%##########################################################################
% methods:
%--------------------------------------------------------------------------
% *) function obj = KiLCA_backgrounddata(path)
% *) function plotB(obj, type, varargin)
%##########################################################################

%author:   Philipp Ulbl
%created:  19.02.2019
%modified: 22.08.2019
        
    properties
        path    %path of background-data
    end
    
    properties (Dependent)
        R       %radius
        b0      %absolute magnetic field b0 = sqrt(b0th^2+b0z^2)
        b0th    %magnetic field in theta direction
        b0z     %magnetic field in Z direction
        dPHI0   %
        Er_i    %radial electric field (from input profiles)
        hth     %
        hz      %
        j0p_m   %
        j0s_m   %
        j0th    %
        j0th_m  %
        j0z     %
        j0z_m   %
        ne_m    % 
        ne_p    %
        n_i     %number density (from input profiles)
        ni_m    %
        ni_p    %
        nue     %
        nui     %
        q_i     %safety factor (from input profiles)
        Te_i    %electron temobjature (from input profiles)
        Te_m    %
        Ti_i    %ion temobjature (from input profiles)
        Ti_m    %
        Vpe_m   %
        Vpe_p   %
        Vpi_m   %
        Vpi_p   %
        Vs0e_m  %
        Vs0i_m  %
        Vse_m   %
        Vsi_m   %
        Vte_p   %
        Vth_i   %theta velocity (from input profiles)
        Vth_m   %
        Vti_p   %
        Vz_i    %Z velocity (from input profiles)
        Vz_m    %
    end
    
    %get Access for dependent proobjties
    methods
       function q = get.R(obj)
           q = get_background_quantity(obj, 'b0');
           q = q(:, 1);
       end
       function q = get.b0(obj)
           q = get_background_quantity(obj, 'b0');
       end
       function q = get.b0th(obj)
           q = get_background_quantity(obj, 'b0th');
       end
       function q = get.b0z(obj)
           q = get_background_quantity(obj, 'b0z');
       end
       function q = get.dPHI0(obj)
           q = get_background_quantity(obj, 'dPHI0');
       end
       function q = get.Er_i(obj)
           q = get_background_quantity(obj, 'Er_i');
       end
       function q = get.hth(obj)
           q = get_background_quantity(obj, 'hth');
       end
       function q = get.hz(obj)
           q = get_background_quantity(obj, 'hz');
       end
       function q = get.j0p_m(obj)
           q = get_background_quantity(obj, 'j0p_m');
       end
       function q = get.j0s_m(obj)
           q = get_background_quantity(obj, 'j0s_m');
       end
       function q = get.j0th(obj)
           q = get_background_quantity(obj, 'j0th');
       end
       function q = get.j0th_m(obj)
           q = get_background_quantity(obj, 'j0th_m');
       end
       function q = get.j0z(obj)
           q = get_background_quantity(obj, 'j0z');
       end
       function q = get.j0z_m(obj)
           q = get_background_quantity(obj, 'j0z_m');
       end
       function q = get.ne_m(obj)
           q = get_background_quantity(obj, 'ne_m');
       end
       function q = get.ne_p(obj)
           q = get_background_quantity(obj, 'ne_p');
       end
       function q = get.n_i(obj)
           q = get_background_quantity(obj, 'n_i');
       end
       function q = get.ni_m(obj)
           q = get_background_quantity(obj, 'ni_m');
       end
       function q = get.ni_p(obj)
           q = get_background_quantity(obj, 'ni_p');
       end
       function q = get.nue(obj)
           q = get_background_quantity(obj, 'nue');
       end
       function q = get.nui(obj)
           q = get_background_quantity(obj, 'nui');
       end
       function q = get.q_i(obj)
           q = get_background_quantity(obj, 'q_i');
       end
       function q = get.Te_i(obj)
           q = get_background_quantity(obj, 'Te_i');
       end
       function q = get.Te_m(obj)
           q = get_background_quantity(obj, 'Te_m');
       end
       function q = get.Ti_i(obj)
           q = get_background_quantity(obj, 'Ti_i');
       end
       function q = get.Ti_m(obj)
           q = get_background_quantity(obj, 'Ti_m');
       end
       function q = get.Vpe_m(obj)
           q = get_background_quantity(obj, 'Vpe_m');
       end
       function q = get.Vpe_p(obj)
           q = get_background_quantity(obj, 'Vpe_p');
       end
       function q = get.Vpi_m(obj)
           q = get_background_quantity(obj, 'Vpi_m');
       end
       function q = get.Vpi_p(obj)
           q = get_background_quantity(obj, 'Vpi_p');
       end
       function q = get.Vs0e_m(obj)
           q = get_background_quantity(obj, 'Vs0e_m');
       end
       function q = get.Vs0i_m(obj)
           q = get_background_quantity(obj, 'Vs0i_m');
       end
       function q = get.Vse_m(obj)
           q = get_background_quantity(obj, 'Vse_m');
       end
       function q = get.Vsi_m(obj)
           q = get_background_quantity(obj, 'Vsi_m');
       end
       function q = get.Vte_p(obj)
           q = get_background_quantity(obj, 'Vte_p');
       end
       function q = get.Vth_i(obj)
           q = get_background_quantity(obj, 'Vth_i');
       end
       function q = get.Vth_m(obj)
           q = get_background_quantity(obj, 'Vth_m');
       end
       function q = get.Vti_p(obj)
           q = get_background_quantity(obj, 'Vti_p');
       end
       function q = get.Vz_i(obj)
           q = get_background_quantity(obj, 'Vz_i');
       end
       function q = get.Vz_m(obj)
           q = get_background_quantity(obj, 'Vz_m');
       end
    end
    
    methods (Access = public)
        function obj = KiLCA_data_background(path)
            %##############################################################
            %function obj = KiLCA_data_background(path)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % constructor of the class. needs the path of the
            % background-data folder.
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % path  ... path of background-data folder
            %##############################################################
              
            obj.path = path;
        end
        
        function plotB(obj, varargin)
            %##############################################################
            %function plotB(obj, varargin)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % plots magnetic field components as 3 subplots in a row
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % varargin  ... plot arguments
            %##############################################################
            
            figure('units', 'normalized', 'outerposition', [0, 0, 1, 1]);
            %call superclass method
            plot_triple(obj, 'b0', 'b0th', 'b0z', 'B / G', 'Re', varargin{:});
        end
    end
    
    methods (Access = private)
        function q = get_background_quantity(obj, file)
            %returns the quantity of the specified file as a matrix with 
            %the radius in the 1st column and the quantity in the 2nd
            q = importdata([obj.path, file, '.dat']);
        end
    end
end