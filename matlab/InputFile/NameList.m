classdef NameList < dynamicprops
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

    methods
        function create(obj, raw)
            %##############################################################
            %function create(obj, raw)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % Creates namelist object out of cell array with lines from the
            % namelistfile in each cell.
            %##############################################################
            % raw   ... cell array with 1 line of namelist file per cell
            %##############################################################
            
            %go through raw and add properties to class
            for k = 1:numel(raw)
                %split by =
                sp = strsplit(raw{k}, '=');
                %get property name
                prop = strrep(strtrim(sp{1}), ' ', '_'); %replace spaces
                %get value as string
                val = strtrim(strrep(sp{2}, ',', '')); %remove comma
                
                %remove comment on the right and trim
                val = strsplit(val, '!');
                val = strtrim(val{1});
                
                %replace t or f by true or false if len=1 (string has len
                %minimum 3 because of "")
                if(numel(val) == 1)
                    val = strrep(val, 't', '.true.');
                    val = strrep(val, 'f', '.false.');
                end
                
                %ADD ARRAY IMPLEMENTATION
                
                %case 1: string
                if(contains(val, '''') || contains(val, '"'))
                    %remove ', "
                    val = strrep(val, '''', '');
                    val = strrep(val, '"', '');
                    
                %case 2: boolean
                elseif(strcmp(val(1), '.') && strcmp(val(end), '.'))
                    val = logical(strcmp(val, '.true.'));
                    
                %case 3: number
                else
                    %convert to number (replace d by e if exist)
                    val = sscanf(strrep(val, 'd', 'e'), '%g');
                end
                addprop(obj, prop);
                obj.(prop) = val;
            end
        end
        
        function write(obj, fname)
            %##############################################################
            %function write(obj, fname)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % Writes namelist object to file. Appends to existing entries
            % in file. Does not write the start and end symbol of the
            % namelist.
            %##############################################################
            % fname  ... name of file to write object to
            %##############################################################
            
            %open file to append
            fid = fopen(fname, 'a');
            %write all properties
            propnames = properties(obj);
            for k = 1:numel(propnames)
                
                %q is the property value
                q = obj.(propnames{k});
                
                %case boolean
                if(isa(q, 'logical'))
                    l = {'.false.', '.true.'};
                    s = l{q+1};
                    
                %case numbers
                elseif(isa(q, 'double'))
                    s = sprintf('%g', q);
                    
                %case strings
                else
                    s = ['''', q, ''''];
                end
                
                %write string to file
                fprintf(fid, '%s\n', [propnames{k}, ' = ', s, ' ,']);
            end
            %close file
            fclose(fid);
        end
    end
end