classdef (Abstract) hdf5_output < handle
%##########################################################################
%classdef (Abstract) hdf5_output < handle
%##########################################################################
% description of class:
%--------------------------------------------------------------------------
% This class is a prototype for hdf5 output. It forces to give a
% description and a physical unit when saving a quantity to a hdf5 file.
% If data is complex, real and imag part are split and concatinated along
% the outermost dimension of the array.
%##########################################################################
% properties:
%--------------------------------------------------------------------------
%##########################################################################
% methods:
%--------------------------------------------------------------------------
% *) writeHDF5(obj, fname, loc, quant, desc, unit)
%########################################################################## 

%author:   Philipp Ulbl
%created:  xx.03.2020

    properties(Access = private)
        ENC = 'UTF8';
    end

    methods(Access = protected)
        
        function writeHDF5(obj, fname, loc, quant, desc, unit)
            %##############################################################
            %function writeHDF5(obj, fname, loc, quant, desc, unit)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % writes a single entry to hdf5 file.
            % If data is complex, real and imag part are split and 
            % concatinated along the outermost dimension of the array.
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % fname ... filename
            % loc   ... location in tree
            % quant ... name of quantitiy (class property)
            % desc  ... short description of quantity
            % unit  ... physical unit of quantity
            %############################################################## 
            
            %name of entry
            name = [loc, quant];
            %data to store
            data = obj.(quant);
            
            %check for string and modify data if yes
            [data, flag_string] = obj.checkString(data);
            
            %check for complex and modify data if yes
            [data, flag_complex, dim] = obj.checkComplex(data);
            
            %write entry
            obj.writeHDF5Entry(fname, name, data, desc, unit);
            
            %if type is string add attribute
            if(flag_string == true)
                h5writeatt(fname, name, 'stringvalued', ['Encoding: ', obj.ENC]);
            end
            %if type is complex add attribute
            if(flag_complex == true)
                h5writeatt(fname, name, 'complexvalued', ['along dim = ', num2str(dim)]);
            end
        end
    end
    
    methods(Access = private)
        
        function [data, flag] = checkString(obj, data)
            
            %if data is string or char convert to utf8
            if(ischar(data) || isstring(data))
                
                data = unicode2native(data, obj.ENC);
                
                flag = true;
            else
                flag = false;
            end
        end
        
        function [data, flag, dim] = checkComplex(obj, data)
            
            %if data is complex, split real and imag
            if(any(imag(data) ~= 0, 'all'))
                %for row vector cat along inner dimension
                if(numel(size(data)) == 2 && size(data, 1) == 1)
                    dim = 1;
                %otherwise along outermost dim
                else
                    dim = sum(size(data) > 1) + 1;
                end
                %add imag to outermost dimension
                data = cat(dim, real(data), imag(data));
                flag = true;
            else
                dim = [];
                flag = false;
            end
        end
        
        function writeHDF5Entry(obj, fname, name, data, desc, unit)
            %create entry with datatype of quantity
            h5create(fname, name, size(data), 'Datatype', class(data));
            %write quantity
            h5write(fname, name, data);
            %write description as attribute
            h5writeatt(fname, name, 'decription', desc);
            %write unit as attribute
            h5writeatt(fname, name, 'unit', unit);
        end
    end
end