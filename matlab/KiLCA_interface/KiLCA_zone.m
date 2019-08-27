classdef KiLCA_zone < KiLCA_prototype_input
%classdef KiLCA_zone
%##########################################################################
% description of class:
%--------------------------------------------------------------------------
% class containing the information of a zone in KiLCA. Corresponds to a
% zoneX.in file for KiLCA.
%##########################################################################
% properties:
%--------------------------------------------------------------------------
% *) r1, typeBC1, model, modelvers, typeBC2, r2
% *) vacuum, imhd, flre
% READONLY:
% *) INDICES, BLUEPRINT, READY
%##########################################################################
% methods:
%--------------------------------------------------------------------------
% *) function obj = KiLCA_zone(r1, b1, m, r2, b2)
% *) function c = plain(obj)
% *) function i = indices(obj)
%##########################################################################

%author:   Philipp Ulbl
%created:  05.02.2019
%modified: 21.08.2019
    
    properties (SetAccess = 'private')
        INDICES         %indices of parameters in blueprint files
        BLUEPRINT       %name of blueprint file
        READY = false;  %flag: ready to run
    end
    
    properties
        r1 = [];        %r1 - minimum radius of the zone (plasma radius)
        typeBC1 = [];   %type of BC at r1 (center, infinity, interface, antenna)
        model =[];      %type of the plasma model (vacuum, medium, imhd, rmhd, flre)
        modelvers = 0;  %code version for the model: MHD model (0 - incompressible and flowless, 1 - compressible with flows)
        typeBC2 = [];   %type of BC at r2 (center, infinity, interface, antenna)
        r2 = [];        %r2 - maximum radius of the zone (first wall boundary)
        
        vacuum          %contains vacuum information about the zone
        imhd            %contains imhd information about the zone
        flre            %contains flre information about the zone
    end
    
    methods
        function obj = KiLCA_zone(r1, b1, m, r2, b2)
            %##############################################################
            %function obj = KiLCA_zone(r1, b1, m, r2, b2)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % constructor of the class. needs the specification of the zone
            % (borders and model) to work.
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % r1    ... inner radius 
            % b1    ... inner boundary type
            % m     ... model between borders
            % r2    ... outer radius
            % b2    ... outer boundary type
            %##############################################################
            
            %check input
            if ~isnumeric(r1) || ~isnumeric(r2)
                error('r1, r2 must be numeric.')
            end
            if ~ischar('b1') || ~ischar('b2') || ~ischar('m')
                error('b1, b2, m must be char.')
            end
            
            obj.r1 = r1;
            obj.typeBC1 = b1;
            obj.model = m;
            obj.r2 = r2;
            obj.typeBC2 = b2;
            
%             if strcmp(m, 'imhd') || strcmp(m, 'rmhd')
%                error('imhd or rmhd support not yet available.') 
%             end
                        
            %choose which class represents the medium
            if strcmp(m, 'vacuum')
            	obj.vacuum = KiLCA_zone_vacuum();
                o = obj.vacuum;

%                 THIS WAS USED SOMEDAY BEFORE TO SIMULATE RESISTIVE WALL
%                 MAY BE OF USE IN FUTURE.
%                 %set idealwall or medium conductivity
%                 if strcmp(b2, 'idealwall') || strcmp(m, 'medium')
%                    obj.vacuum.sigma(1) = 1.3e16;
%                 end

            elseif strcmp(m, 'imhd')
                obj.imhd = KiLCA_zone_imhd();
                o = obj.imhd;
            elseif strcmp(m, 'flre')
                obj.flre = KiLCA_zone_flre();    
                o = obj.flre;
            end
            
            %set INDICES and BLUEPRINT property
            obj.INDICES   = o.INDICES;
            obj.BLUEPRINT = o.BLUEPRINT;
            
            obj.READY = true;
        end
        
        function c = plain(obj)
            %##############################################################
            %function c = plain(obj)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % returns the class as an cell-array with the name of the
            % properties in each cell. properties which are structs will be
            % split into their elements.
            % *extends superclass .plain() method
            %##############################################################
            % output:
            %--------------------------------------------------------------
            % c     ... cell array
            %##############################################################
            
            %use superclass method
            c1 = plain@KiLCA_prototype_input(obj);
            c1(end) = []; % delete flre class
            c1(end) = []; % delete imhd class
            c1(end) = []; % delete vacuum class
            
            %extend array by cell array of the used zone class
            if strcmp(obj.model, 'vacuum') 
                c2 = obj.vacuum.plain();
            elseif strcmp(obj.model, 'imhd')
                c2 = obj.imhd.plain();
            elseif strcmp(obj.model, 'flre')
                c2 = obj.flre.plain();
            else
                error('type of zone not supported (only vacuum, imhd, flre).')
            end
            
            c = [c1, c2];
        end
    end
end