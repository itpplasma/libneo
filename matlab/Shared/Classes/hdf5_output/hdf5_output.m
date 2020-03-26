classdef (Abstract) hdf5_output < handle
%##########################################################################
%classdef (Abstract) hdf5_output < handle
%##########################################################################
% description of class:
%--------------------------------------------------------------------------
% This class is a prototype for hdf5 output. It forces to give a
% description and a physical unit when saving a quantity to a hdf5 file.
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
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % fname ... filename
            % loc   ... location in tree
            % quant ... name of quantitiy (class property)
            % desc  ... short description of quantity
            % unit  ... physical unit of quantity
            %############################################################## 
            
            %create entry with datatype of quantity
            h5create(fname, [loc, quant], size(obj.(quant)), 'Datatype', class(obj.(quant)));
            %write quantity
            h5write(fname, [loc, quant], obj.(quant));
            %write description as attribute
            h5writeatt(fname, [loc, quant], 'decription', desc);
            %write unit as attribute
            h5writeatt(fname, [loc, quant], 'unit', unit);
        end
    end
end