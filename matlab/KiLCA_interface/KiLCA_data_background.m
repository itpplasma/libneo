classdef KiLCA_data_background < KiLCA_prototype_output & hdf5_output
%classdef KiLCA_data_background < KiLCA_prototype_output & hdf5_output
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
% *) j0p_m, j0s_m, j0th, j0th_m in units of c=1 !!!
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
% *) function plotProf(obj, modes, varargin)
%##########################################################################

%author:   Philipp Ulbl
%created:  19.02.2019
%modified: 04.11.2019
        
    properties
        path    %path of background-data
    end
    
    properties (Access = private)
        c = 2.99792458e10;
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
        Te_i    %electron temperature (from input profiles)
        Te_m    %
        Ti_i    %ion temperature (from input profiles)
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
        
        %self defined
        j0
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
           q(:, 2) = q(:, 2) ./ obj.c;
       end
       function q = get.j0s_m(obj)
           q = get_background_quantity(obj, 'j0s_m');
           q(:, 2) = q(:, 2) ./ obj.c;
       end
       function q = get.j0th(obj)
           q = get_background_quantity(obj, 'j0th');
           q(:, 2) = q(:, 2) ./ obj.c;
       end
       function q = get.j0th_m(obj)
           q = get_background_quantity(obj, 'j0th_m');
           q(:, 2) = q(:, 2) ./ obj.c;
       end
       function q = get.j0z(obj)
           q = get_background_quantity(obj, 'j0z');
           q(:, 2) = q(:, 2) ./ obj.c;
       end
       function q = get.j0z_m(obj)
           q = get_background_quantity(obj, 'j0z_m');
           q(:, 2) = q(:, 2) ./ obj.c;
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
       
       %self defined
       function q = get.j0(obj)
           q = hypot(obj.j0th, obj.j0z);
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
        
        function plotJ(obj, varargin)
            %##############################################################
            %function plotJ(obj, varargin)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % plots current density components as 3 subplots in a row
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % varargin  ... plot arguments
            %##############################################################
            
            figure('units', 'normalized', 'outerposition', [0, 0, 1, 1]);
            %call superclass method
            plot_triple(obj, 'j0', 'j0th_m', 'j0z_m', 'J / statA cm^{-2}', 'Re', varargin{:});
        end
        
        function plotProf(obj, modes, varargin)
            %##############################################################
            %function plotProf(obj, modes, varargin)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % plots all relevant profiles together with the location of the
            % resonant surfaces for all modes given in the input modes.
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % modes     ... KiLCA_modes class
            % varargin  ... plot arguments
            %##############################################################
            
            figure('units', 'normalized', 'outerposition', [0, 0, 1, 1]);
            
            hold on
            
            %plot rescaled perpendicular electron fluid velocity
            vmax = -floor(log10(max(abs(obj.Vs0e_m(:, 2)))));
            vs = sign(obj.Vs0e_m(1, 2));
            if(vs > 0) 
                s = ''; 
            else
                s = '-';
            end
            plot(obj.Vs0e_m(:, 1), vs .* 10^(vmax) .* obj.Vs0e_m(:, 2), ':', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 3, 'DisplayName', [s, 'v_{e,\perp} / 10^{', num2str(-vmax), '} cm s^{-1}']);
            
            %plot positive q
            qs = sign(obj.q_i(1, 2));
            if(qs > 0) 
                s = ''; 
            else
                s = '-';
            end
            plot(obj.q_i(:, 1), qs .* obj.q_i(:, 2), 'k', 'LineWidth', 3, 'DisplayName', [s, 'q']);
            
            %plot rescaled density
            nmax = -floor(log10(max(abs(obj.n_i(:, 2)))));
            plot(obj.n_i(:, 1), 10^(nmax) .* obj.n_i(:, 2), '-.' , 'Color', [0.3010 0.7450 0.9330], 'LineWidth', 3,  'DisplayName', ['n / 10^{', num2str(-nmax), '} cm^{-3}']);
            
            %plot rescaled temperatures
            Tmax = -floor(log10(max(abs(obj.Te_i(:, 2)))));
            plot(obj.Ti_i(:, 1), 10^(Tmax) .* obj.Ti_i(:, 2), '--', 'Color', [0 0.4470 0.7410], 'LineWidth', 3,  'DisplayName', ['T_i / 10^{', num2str(-Tmax), '} keV']);
            plot(obj.Te_i(:, 1), 10^(Tmax) .* obj.Te_i(:, 2), '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 3,  'DisplayName', ['T_e / 10^{', num2str(-Tmax), '} keV']);
            
            %plot lines for modes
            for k=1:numel(modes.m)
                %calculate position of resonant surface
                rres = interp1(obj.q_i(:, 2), obj.q_i(:, 1), -modes.m(k)/modes.n(k));
                %plot line with resonant surface location
                plot([rres, rres], ylim, ':m', 'LineWidth', 2, 'HandleVisibility','off');
                %show modenumber on plot near the line (on pos/neg y alternating)
                text(rres, 0.2, vec2str([modes.m(k),modes.n(k)], '%d'))
            end
            
            %plot zero line
            plot(xlim, [0, 0], '-k', 'LineWidth', 1, 'HandleVisibility','off');
            
            hold off
            
            legend('Location', 'NorthWest')
            
            xlabel('r / cm')
            ylabel('rescaled profiles')
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
            
            obj.writeHDF5(fname, loc, 'R', 'small radius vector', 'cm');
            
            obj.writeHDF5(fname, loc, 'b0', 'absolute magnetic field', 'G');
            obj.writeHDF5(fname, loc, 'b0th', 'poloidal magnetic field', 'G');
            obj.writeHDF5(fname, loc, 'b0z', 'toroidal magnetic field', 'G');
            
            obj.writeHDF5(fname, loc, 'j0', 'absolute current density', 'statA cm^{-2} c=1');
            obj.writeHDF5(fname, loc, 'j0th', 'poloidal current density', 'statA cm^{-2} c=1');
            obj.writeHDF5(fname, loc, 'j0z', 'toroidal current density', 'statA cm^{-2} c=1');
        end
    end
    
    methods (Access = protected)
        function writeHDF5(obj, fname, loc, quant, desc, unit)
            %##############################################################
            %function writeHDF5(obj, fname, loc, quant, desc, unit)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % writes a single entry to hdf5 file. overrides superclass
            % method.
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % fname ... filename
            % loc   ... location in tree
            % quant ... name of quantitiy (class property)
            % desc  ... short description of quantity
            % unit  ... physical unit of quantity
            %############################################################## 
            
            %if quantity is vector call superclass method
            if(sum(size(obj.(quant)) == 1) == 1)
                writeHDF5@hdf5_output(obj, fname, loc, quant, desc, unit);
            %if quantity is matrix of type [R; quant]
            else
                data = obj.(quant);
                %create entry with datatype of quantity
                h5create(fname, [loc, quant], size(data), 'Datatype', class(data));
                %write quantity
                h5write(fname, [loc, quant], data);
                %write description as attribute
                h5writeatt(fname, [loc, quant], 'decription', desc);
                %write unit as attribute
                h5writeatt(fname, [loc, quant], 'unit', unit);
            end
        end
    end
    
    methods (Access = private)
        function q = get_background_quantity(obj, file)
            %returns the quantity of the specified file as a matrix with 
            %the radius in the 1st column and the quantity in the 2nd
            q = load([obj.path, file, '.dat']);
        end
    end
end
