classdef KiLCA_antenna < KiLCA_prototype_input
%classdef KiLCA_antenna
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
%modified: 21.08.2019

    properties (SetAccess = 'private')
        INDICES = 2:8;              %indices of parameters in blueprint files
        BLUEPRINT = 'antenna.in';   %name of blueprint file
        READY = false;              %flag: ready to run
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
    end
end