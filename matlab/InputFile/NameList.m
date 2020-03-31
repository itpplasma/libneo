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
                %metaproperty object of ne property
                metaprop = [];
                %try to parse line
                try
                    %split by =
                    sp = strsplit(raw{k}, '=');
                    %get property name
                    prop = strrep(strtrim(sp{1}), ' ', '_'); %replace spaces
                    %get value as string
                    val = strtrim(strrep(sp{2}, ',', '')); %remove comma

                    %remove comment on the right and trim
                    val = strsplit(val, '!');
                    val = strtrim(val{1});

                    %New Method: 
                    %convert namelist expression to MATLAB expression

                    %get propertyname
                    if(contains(prop, '(') && contains(prop, ')')) %array
                        p = strsplit(prop, '(');
                        %propertyname is the part before the (...)
                        propname = p{1};
                        %array to cell array to allow empty values
                        val = ['{', val, '}'];
                        %create property string for access in MATLAB
                        prop = ['("', strrep(prop, '(', '")(')];
                        %init value for array
                        init = '{}';
                    else %scalars
                        propname = prop;
                        %create property string for access in MATLAB
                        prop = ['("', prop, '")'];
                        %init value for scalar
                        init = '[]';
                    end

                    %handle decimals and booleans
                    if(~contains(val, '''') && ~contains(val, '"'))

                        %replace fortran d by e
                        val = strrep(val, 'd', 'e');

                        %bring all boolean variants in t,f form
                        %we do this because:
                        %if you replace all t by true, true becomes truerue.
                        val = strrep(val, '.FALSE.', 'f');
                        val = strrep(val, '.TRUE.', 't');
                        val = strrep(val, '.False.', 'f');
                        val = strrep(val, '.True.', 't');
                        val = strrep(val, '.false.', 'f');
                        val = strrep(val, '.true.', 't');
                        val = strrep(val, 'false', 'f');
                        val = strrep(val, 'true', 't');

                        %replace t,f for MATLAB expressions
                        val = strrep(val, 't', 'true');
                        val = strrep(val, 'f', 'false');
                    end

                    %add property if not existend and init
                    if(~any(strcmp(properties(obj), propname)))
                        metaprop = addprop(obj, propname);
                        eval(['obj.("', propname, '") = ', init, ';']);
                    end

                    %set property value
                    eval(['obj.', prop, ' = ', val, ';']);
                    
                catch
                    warning(['could not understand line:\n', raw{k}, '\nproperty skipped.'], 'all')
                    %delete property if it has been created
                    if(~isempty(metaprop))
                        delete(metaprop);
                    end
                end
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