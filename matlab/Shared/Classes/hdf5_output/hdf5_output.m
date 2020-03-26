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
            
            %data to store
            data = obj.(quant);
            
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
                flag_complex = true;
            else
                flag_complex = false;
            end
            
            %create entry with datatype of quantity
            h5create(fname, [loc, quant], size(data), 'Datatype', class(data));
            %write quantity
            h5write(fname, [loc, quant], data);
            %write description as attribute
            h5writeatt(fname, [loc, quant], 'decription', desc);
            %write unit as attribute
            h5writeatt(fname, [loc, quant], 'unit', unit);
            
            %if type is complex add attribute
            if(flag_complex == true)
                h5writeatt(fname, [loc, quant], 'complexvalued', ['along dim = ', num2str(dim)]);
            end
        end
    end
end