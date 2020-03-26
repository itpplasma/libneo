classdef (Abstract) blueprint < handle
%classdef (Abstract) blueprint < handle
%##########################################################################
% description of class:
%--------------------------------------------------------------------------
% This class is a prototype for all classes that are interfaces to input
% files that need "Blueprints". This means, they take existing files that
% work (=Blueprints), change the numbers and write the new file.
%##########################################################################
% properties:
%--------------------------------------------------------------------------
% *) READY
%##########################################################################
% methods:
%--------------------------------------------------------------------------
% *) function c = plain(obj)
% *) function write(obj, path_from, path_to, zone_num)
%##########################################################################

%author:   Philipp Ulbl
%created:  07.01.2020
    
    properties (SetAccess = 'protected')
        READY = false;  %flag: ready to run
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
        
        function write(obj, path_from, path_to)
            %##############################################################
            %function write(obj, path_from, path_to)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % writes properties of the class into input files.
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % path_from  ... path where the blueprint is from
            % path_to    ... path where the input file will be written
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
            %change options in cell array. sep = 2xspaces
            raw = change_opts(raw, obj.INDICES, obj.plain(), '  ');
            
            %for not KiLCA zone classes
            save_file(raw, [path_to, obj.BLUEPRINT]);
        end
    end
end
