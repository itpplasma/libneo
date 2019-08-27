classdef (Abstract) KiLCA_prototype_input < handle
%classdef (Abstract) KiLCA_prototype_input < handle
%##########################################################################
% description of class:
%--------------------------------------------------------------------------
% prototype class for all KiLCA input classes.
%##########################################################################
% properties:
%--------------------------------------------------------------------------
%##########################################################################
% methods:
%--------------------------------------------------------------------------
% *) function c = plain(obj)
% *) function write(obj, path_from, path_to, zone_num)
%##########################################################################

%author:   Philipp Ulbl
%created:  21.08.2019
%modified: 21.08.2019
    
    properties
    end
    
    methods (Access = 'public')
        function c = plain(obj)
            %##############################################################
            %function c = plain(obj)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % returns the class as an cell-array with the name of the
            % properties in each cell. properties which are structs will be
            % split into their elements.
            %##############################################################
            % output:
            %--------------------------------------------------------------
            % c     ... cell array
            %##############################################################
            
            c = classprop2cell(obj);
        end
        
        function write(obj, path_from, path_to, zone_num)
            %##############################################################
            %function write(obj, path_from, path_to, zone_num)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % writes properties of the class into a KiLCA input file.
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % path_from  ... path where the blueprint is from
            % path_to    ... path where the input file will be written
            % zone_num   ... number of the zone (needed for KiLCA_zone)
            %##############################################################
            
            %check if class is meant to be written
            if isprop(obj, 'INDICES') && obj.READY == 0
                error('class is not able to be written into a file.')
            end
            
            %check if class is ready to be written
            if isprop(obj, 'READY') && obj.READY == 0
                error('class is not ready to run.')
            end
            
            %read in blueprint
            raw = read_in([path_from, obj.BLUEPRINT]);
            %change options in cell array
            raw = change_opts(raw, obj.INDICES, obj.plain());
            
            %save new file
            if ~isa(obj, 'KiLCA_zone')
                %for not KiLCA zone classes
                save_file(raw, [path_to, obj.BLUEPRINT]);
            else
                %for KiLCA zone classes
                if nargin < 4 || isempty(zone_num)
                    error('zone_num must be given for KiLCA_zone classes.')
                end
                %replace zone in BLUEPRINT name with zone#
                save_file(raw, [path_to, ...
                    strrep(obj.BLUEPRINT, 'zone', ['zone_', num2str(zone_num)])]);
            end
        end
    end
end