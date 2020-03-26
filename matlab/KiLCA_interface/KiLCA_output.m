classdef KiLCA_output < handle & blueprint
%classdef KiLCA_output < handle & blueprint
%##########################################################################
% description of class:
%--------------------------------------------------------------------------
% class containing the information of the output in KiLCA.
%
% General run settings: 0 - to skip, 1 - to calculate only,
% 2 - to calculate and save to a disk
%##########################################################################
% properties: (all have defaults!)
%--------------------------------------------------------------------------
% *) backdata, lindata, varquant, disp
% *) flag_deb
% *) curdenspet, abspow, dispow, kinfluxr, poyntfluxr, totfluxr,
%    numdenspet, lortorq
% READONLY:
% *) INDICES, BLUEPRINT, READY
%##########################################################################
% methods:
%--------------------------------------------------------------------------
% *) function obj = KiLCA_output()
% *) function c = plain(obj)
%##########################################################################
% comments:
%--------------------------------------------------------------------------
% If apropriate the quantity is computed for ions (i), electrons (e) and total (t = i + e).
% 
% If appropriate the quantity is computed for several types of the current, j0 (0) - basic current densit
% 
% Bear in mind that some quantities are dependent on another ones!
%##########################################################################

%author:   Philipp Ulbl
%created:  08.08.2019
%modified: 21.08.2019

    properties (Transient, SetAccess = 'protected')
        INDICES = [2:5, 8, 11:18];  %indices of parameters in blueprint files
        BLUEPRINT = 'output.in';    %name of blueprint file
        SEP = '#'
    end
    
    properties
        backdata    = 2     %background data: default=2
        lindata     = 2     %linear data: default=2
        varquant    = 2     %various quantities (see below): default=2
        disp        = 0     %dispersion: default=0
        flag_deb    = 0     %flag for debugging
        curdenspet  = 1     %flag: current density perturbation
        abspow      = 1     %flag: absorbed power (density and integrated over cylinder volume)
        dispow      = 1     %flag: dissipated power (density and integrated over cylinder volume)
        kinfluxr    = 1     %flag: r-component of kinetic flux out of the cylinder volume
        poyntfluxr  = 1     %flag: r-component of Poynting flux out of the cylinder volume
        totfluxr    = 1     %flag: r-component of total flux out of the cylinder volume
        numdenspet  = 1     %flag: number density perturbation
        lortorq     = 1     %flag: Lorentz torque (density and integrated over cylinder volume)
    end
    
    methods
        function obj = KiLCA_output()
            %##############################################################
            %function obj = KiLCA_output()
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % constructor of the output class. empty.
            %##############################################################
            
            obj.READY = true;
        end
    end
end