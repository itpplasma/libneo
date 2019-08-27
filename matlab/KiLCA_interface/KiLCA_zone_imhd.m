classdef KiLCA_zone_imhd < KiLCA_prototype_input
%classdef KiLCA_zone_imhd
%##########################################################################
% description of class:
%--------------------------------------------------------------------------
% class representing a imhd zone in KiLCA
%##########################################################################
% properties:
%--------------------------------------------------------------------------
% *) sigma
% *) rgrid_maxdim, relacc, absacc
% *) polydeg, sparse_relacc, sparse_absacc, maxgridstep
% *) flag_deb
% READONLY:
% *) INDICES, BLUEPRINT, READY
%##########################################################################
% methods:
%--------------------------------------------------------------------------
% *) function obj = KiLCA_zone_imhd()
% *) function c = plain(obj)
%##########################################################################
    
    %author:   Philipp Ulbl
    %created:  21.08.2019
    %modified: 21.08.2019
    

    properties (SetAccess = 'private')
        INDICES   = [2:7, 10:12, 15:18, 21]; %indices of parameters in blueprint files
        BLUEPRINT = 'zone_imhd.in';            %name of blueprint file
    end
    
    properties
        %ME solver settings:
        
        rgrid_maxdim    = 1e6    %max dimension of the radial grid for the solution: default=1e5
        relacc          = 1e-12  %relative accuracy of the solution: default=1e-8
        absacc          = 1e-12  %absolute accuracy of the solution: default=1e-8
        
        %ME solution space out settings:
        
        polydeg         = 3      %degree of the polynomial used to space out the solution (by checking the accuracy
        sparse_relacc   = 1e-8   %relative accuracy of the sparse solution: default=1e-8
        sparse_absacc   = 1e-8   %absolute accuracy of the sparse solution: default=1e-8
        maxgridstep     = 0.1    %max grid step in the solution: default=0.1
        
        flag_deb        = 0      %flag for debugging mode (additional checks are performed in the code): default=1
    end
    
    methods
        function obj = KiLCA_zone_imhd()
            %##############################################################
            %function obj = KiLCA_zone_imhd()
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % constructor of the vacuum zone class. empty.
            %##############################################################
        end
    end
end