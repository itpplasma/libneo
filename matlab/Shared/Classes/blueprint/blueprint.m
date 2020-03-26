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
    
    properties (Abstract, Transient, SetAccess = 'protected')
        INDICES
        BLUEPRINT
        SEP
    end
    properties (Transient, SetAccess = 'protected')
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
            
            %check for correct implementation of class
            obj.validate();
            
            %check if class is ready to be written
            if obj.READY == 0
                error('class is not ready to run.')
            end
            
            %read in blueprint
            raw = read_in([path_from, obj.BLUEPRINT]);
            %change options in cell array. sep = 2xspaces
            raw = change_opts(raw, obj.INDICES, obj.plain(), obj.SEP);
            
            %save file
            save_file(raw, [path_to, obj.BLUEPRINT]);
        end
    end
    
    methods (Access = 'public')
        
        function validate(obj)
            %check for correct implementation of class
            
            %get metaclass object of object
            prop = metaclass(obj);
            %get metaclass object of blueprint
            prop_blue = ?blueprint;
            
            %go through all properties in blueprint
            for k = 1:numel(prop_blue.PropertyList)
                %check only for transient
                if (prop_blue.PropertyList(k).Transient == true)
                    %get index of property that coincides with property k
                    ind = arrayfun(@(a) strcmp(a.Name, prop_blue.PropertyList(k).Name), prop.PropertyList);
                    %check if it is transient and throw error if not
                    if any(prop.PropertyList(ind).Transient == false)
                        error('implementations of transient properties of blueprint class must be transient.')
                    end
                end
            end
        end
    end
end