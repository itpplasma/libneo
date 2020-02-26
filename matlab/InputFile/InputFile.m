classdef InputFile < dynamicprops
%classdef NameList < dynamicprops
%##########################################################################
% description:
%--------------------------------------------------------------------------
% This class is used to store a namelist as an object. Properties are
% created automatically within create and can be accessed by their name.
% When writing, the property name and its value will written to the file.
%##########################################################################
% methods:
%--------------------------------------------------------------------------
% function create(obj, raw)
% function write(obj, fname)
%##########################################################################

    
%author: Philipp Ulbl
%created: 25.02.2020

    properties (SetAccess = private)
        path
    end
    
    methods
        function obj = InputFile(path)
            %##############################################################
            %function obj = InputFile(path)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % Constructor of the class. Checks for file existence.
            %##############################################################
            % path  ... path of input file with name
            %##############################################################
            
            %check if file exists
            if(exist(path, 'file') ~= 2)
               error(['file not found in ', path]); 
            end
            
            obj.path = path;
        end
        
        function read(obj)
            %##############################################################
            %function read(obj)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % Read the Input file. For each namelist inside, a property
            % will be created automatically. These properties are of type
            % NameList and contain the namelist properties which can be
            % accessed by their name.
            %##############################################################
            
            %get data
            raw = importdata(obj.path); %CAN BE IMPROVED
            %trim all lines
            raw = cellfun(@(x) strtrim(x), raw, 'UniformOutput', false);
            %delete comments
            raw(cellfun(@(x) strcmp(x(1), '!'), raw, 'UniformOutput', true)) = [];
            
            %get namelists beginnings and ends
            nml1 = cellfun(@(x) strcmp(x(1), '&'), raw, 'UniformOutput', true);
            nml2 = cellfun(@(x) strcmp(x(1), '/'), raw, 'UniformOutput', true);
            
            %check for error
            if(sum(nml1) ~= sum(nml2))
                error('error reading namelist inside file. number of & and / not the same.')
            end
            
            %start and end indices of all namelists in raw
            ind1 = find(nml1, 1, 'first');
            ind2 = find(nml2, 1, 'first');
            
            %add all namelists as single properties
            for k = 1:numel(ind1)
                %get name of namelist without &
                name = raw{ind1(k)}(2:end);
                %add namelist property to class
                addprop(obj, name);
                %create namelist object
                obj.(name) = NameList();
                obj.(name).create(raw(ind1(k)+1:ind2(k)-1)); %without & and / lines
            end
        end
        
        function write(obj, pathto)
            %##############################################################
            %function write(obj, pathto)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % Writes the input file into the specified path + name.
            % NameList properties will be written. Directory structure will
            % be created if not existent. Existent files will be
            % overwritten.
            %##############################################################
            % pathto  ... path of output file with name
            %##############################################################
            
            %get only dir structure
            [dirto,~,~] = fileparts(pathto);
            %create path if not exist
            system(['mkdir -p ', dirto]);
            %remove file if exist (ignore warnings)
            system(['rm ', pathto, ' 2>/dev/null']);
            %iterate through all properties that are not path and write
            propnames = properties(obj);
            for k = 1:numel(propnames)
                
                %skip property if not of type NameList
                if(~isa(obj.(propnames{k}), 'NameList'))
                    continue;
                end
                
                %write beginning of namelist
                fid = fopen(pathto, 'a');
                fprintf(fid, '%s\n', ['&', propnames{k}]);
                fclose(fid);
                
                %write content
                obj.(propnames{k}).write(pathto);
                
                %write end of namelist
                fid = fopen(pathto, 'a');
                fprintf(fid, '%s\n', '/');
                fclose(fid);
            end
        end
    end
end