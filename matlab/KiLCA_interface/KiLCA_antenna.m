classdef KiLCA_antenna < handle & blueprint & hdf5_output
%classdef KiLCA_antenna < handle & blueprint & hdf5_output
%##########################################################################
% description of class:
%--------------------------------------------------------------------------
% class containing the information of the antenna in KiLCA.
%##########################################################################
% properties:
%--------------------------------------------------------------------------
% *) ra, width, I0, freq, nmod, fdeb, feig
% READONLY:
% *) INDICES, BLUEPRINT, READY
%##########################################################################
% methods:
%--------------------------------------------------------------------------
% *) function obj = antenna(r, n)
% *) function c = plain(obj)
%##########################################################################

%author:   Philipp Ulbl
%created:  08.08.2019
%modified: 26.03.2020

    properties (Transient, SetAccess = 'protected')
        INDICES = 2:8;              %indices of parameters in blueprint files
        BLUEPRINT = 'antenna.in';   %name of blueprint file
        SEP = '#'
    end

    properties
        ra                  %small radius (cm) of antenna location (must match to zones.in)
        width = 0.0         %current density layer width; if == 0 then use delta function
        I0 = 4.5e12         %current in the coils (statamp): default=4.5e12 (1.5 kA)
        flab = [1e3, 0]     %complex frequency (1/s) in the laboratory frame: default=(1.0e3, 0.0e0)
        nmod                %number of antenna modes to be used from modes.in file (stored as: (m - poloidal, n - toroidal))
        flag_deb = 0        %flag for debugging (additional checks are performed in the code): default=1
        flag_eig = 0        %flag to solve an eigenmode problem: 0 - no, 1 - yes
    end
    
    methods (Access = 'public')
        function obj = KiLCA_antenna(r, n)
            %##############################################################
            %function obj = KiLCA_antenna(r, n)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % constructor of the antenna class.
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % r     ... antenna position in cm (minor radius)
            % n     ... number of antenna modes
            %##############################################################
            obj.ra = r;
            obj.nmod = n;
            
            obj.READY = true;
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
            
            obj.writeHDF5(fname, loc, 'ra', 'small radius of antenna location', 'cm');
            obj.writeHDF5(fname, loc, 'I0', 'antenna coil current constant', 'statA');
            obj.writeHDF5(fname, loc, 'flab', 'complex frequency in the laboratory frame', 's^{-1}');
        end
    end
end