classdef KiLCA_zone < handle & blueprint & hdf5_output
%classdef KiLCA_zone < handle & blueprint & hdf5_output
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
%modified: 26.03.2020
    
    properties (Transient, SetAccess = 'protected')
        INDICES         %indices of parameters in blueprint files
        BLUEPRINT       %name of blueprint file
        SEP = '#'
    end
    properties (Transient, SetAccess = 'private')
        number          %zone number
    end
    
    properties
        r1 = [];        %r1 - minimum radius of the zone (plasma radius)
        typeBC1 = [];   %type of BC at r1 (center, infinity, interface, antenna, idealwall)
        model =[];      %type of the plasma model (vacuum, medium, imhd, rmhd, flre)
        modelvers = 0;  %code version for the model: MHD model (0 - incompressible and flowless, 1 - compressible with flows)
        typeBC2 = [];   %type of BC at r2 (center, infinity, interface, antenna, idealwall)
        r2 = [];        %r2 - maximum radius of the zone (first wall boundary)
        
        vacuum          %contains vacuum information about the zone
        imhd            %contains imhd information about the zone
        flre            %contains flre information about the zone
    end
    
    properties (Dependent, Transient)
        typeBC1_num;    %numeric type of BC at r1 (center=0, infinity=1, interface=2, antenna=3, idealwall=4)
        model_num;      %numeric type of the plasma model (vacuum=0, medium=1, imhd=2, rmhd=3, flre=4)
        typeBC2_num;    %numeric type of BC at r1 (center=0, infinity=1, interface=2, antenna=3, idealwall=4)
    end
    
    methods %getter
       function q = get.typeBC1_num(obj)
           q = obj.BC_num('typeBC1');
       end
       function q = get.typeBC2_num(obj)
           q = obj.BC_num('typeBC2');
       end
       function q = get.model_num(obj)
           q = obj.MD_num('model');
       end
    end
    methods
        function obj = KiLCA_zone(num, r1, b1, m, r2, b2)
            %##############################################################
            %function obj = KiLCA_zone(num, r1, b1, m, r2, b2)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % constructor of the class. needs the specification of the zone
            % (borders and model) to work.
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % num   ... zone number
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
            
            if strcmp(m, 'imhd') || strcmp(m, 'rmhd')
               error('imhd or rmhd support not yet available.') 
            end
                        
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
            
            obj.number = num;
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
            c1 = plain@blueprint(obj);
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
        
        function write(obj, path_from, path_to)
            %##############################################################
            %function write(obj, path_from, path_to)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % writes properties of the class into input files. Overrides
            % superclass method because output file name is different from
            % input file name.
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % path_from  ... path where the blueprint is from
            % path_to    ... path where the input file will be written
            %##############################################################
            
            %check if class is ready to be written
            if obj.READY == 0
                error('class is not ready to run.')
            end
            
            %read in blueprint
            raw = read_in([path_from, obj.BLUEPRINT]);
            %change options in cell array. sep = 2xspaces
            raw = change_opts(raw, obj.INDICES, obj.plain(), obj.SEP);
            
            %save file
            save_file(raw, [path_to, strrep(obj.BLUEPRINT, 'zone', ['zone_', num2str(obj.number)])]);
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
            
            obj.writeHDF5(fname, loc, 'r1', 'minimum radius of the zone', 'cm');
            obj.writeHDF5(fname, loc, 'r2', 'maximum radius of the zone', 'cm');
            obj.writeHDF5(fname, loc, 'typeBC1_num', 'numeric type of BC at r1', '(center=0, infinity=1, interface=2, antenna=3, idealwall=4)');
            obj.writeHDF5(fname, loc, 'typeBC2_num', 'numeric type of BC at r2', '(center=0, infinity=1, interface=2, antenna=3, idealwall=4)');
            obj.writeHDF5(fname, loc, 'model_num', 'numeric type of the plasma model', '(vacuum=0, medium=1, imhd=2, rmhd=3, flre=4)');
        end
    end
    
    methods (Access = private)
        function q = BC_num(obj, prop)
           %returns numeric type of BC property
           if(strcmp(obj.(prop), 'center'))
               q=0;
           elseif (strcmp(obj.(prop), 'infinity'))
               q=1;
           elseif (strcmp(obj.(prop), 'interface'))
               q=2;
           elseif (strcmp(obj.(prop), 'antenna'))
               q=3;
           elseif (strcmp(obj.(prop), 'idealwall'))
               q=4;
           end
        end
        function q = MD_num(obj, prop)
            %returns numeric type of model property
            if(strcmp(obj.(prop), 'vacuum'))
                q=0;
            elseif (strcmp(obj.(prop), 'medium'))
                q=1;
            elseif (strcmp(obj.(prop), 'imhd'))
                q=2;
            elseif (strcmp(obj.(prop), 'rmhd'))
                q=3;
            elseif (strcmp(obj.(prop), 'flre'))
                q=4;
            end
        end
    end
end