classdef mnDAT < handle
%classdef mnDAT < handle
%##########################################################################
% description of class:
%--------------------------------------------------------------------------
% This class is intended to be used to create mndat files for the use in
% e.g. Patricks MHD code or in codes of Christopher. This file-type 
% represents basically fourier mode profiles of some quantity q, e.g. Br.
% The profiles can be given on a flux surface function or on a radius. The
% latter is basically the representation in KiLCA, thus such files are
% marked as _KiLCA. 
%
% Planned:
% +) add reading of existent mnDAT files
%##########################################################################
% properties (READONLY):
%--------------------------------------------------------------------------
% fname, kilca, n, s, q, mn
%##########################################################################
% methods:
%--------------------------------------------------------------------------
% *) function obj = mnDAT(path, q, i, n, kilca)
% *) function set(obj, s, q, mn)
% *) function write(obj)
%##########################################################################

%author: Philipp Ulbl
%created: 04.11.2019

    
    properties (SetAccess = private)
        FORMAT = ' %22.15e'; %Format string used to write data: CONSTANT!
        
        fname = {}; %name of file with path
        
        kilca;      %boolean: KiLCA type or not
        
        n           %toroidal modenumber associated to this file
        
        s           %flux surface label: psi or r for KiLCA
        q           %safety factor
        mn          %mn data
    end
    
    methods
        function obj = mnDAT(path, I, i, n, kilca)
            %##############################################################
            %function obj = mnDAT(path, I, i, n, kilca)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % creates an instance of the class.
            %##############################################################
            % path  ... path of the file
            % I     ... quantitiy contained in mn. e.g.: B
            % i     ... index of I. e.g.: r (-> B_r)
            % n     ... toroidal modenumber associated to this file
            % kilca ... boolean to indicate if this file is KiLCA-like
            %##############################################################
            
            %splits the extenson from the filename and path if there
            [p, ~, ~] = fileparts(path);
            
            obj.n = n;
            obj.kilca = kilca;
            
            %add extension to file to identify KiLCA-likeness
            if(kilca)
                k = '_KiLCA';
            else
                k = '';
            end
            
            %construct path+name, add extension .dat
            obj.fname = [p, '/', I, 'm', num2str(n), '_', i, k, '.dat'];
            
        end
        
        function set(obj, s, q, mn)
            %##############################################################
            %function set(obj, s, q, mn)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % fills the class with all important properties.
            %##############################################################
            % s     ... flux surface label: psi or r for KiLCA
            % q     ... safety factor
            % mn    ... mn data matrix
            %##############################################################
            
            mndim = size(mn);
            
            if(any(size(s)<1))
                error('s must be at least scalar.');
            end
            
            if(size(q) ~= size(s))
                error('size of q and s must be the same.');
            end
            if(mndim(1) ~= numel(s))
                error('number of rows in mn must be equal to numel of s.');
            end
            if(mod(mndim(2), 2) ~= 0)
                error('mn must be even in columns.');
            end
            
            obj.s  = s;
            obj.q  = q;
            obj.mn = mn;
        end
        
        
        function write(obj)
            %##############################################################
            %function write(obj)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % writes the data in the class into a mnDAT file.
            %##############################################################
            
            %check if path exists
            if ~isfolder(fileparts(obj.fname))
                mkdir(fileparts(obj.fname));
            end
            
            %check if file exists and delete
            if isfile(obj.fname)
                delete(obj.fname);
                warning(['existing file <', obj.fname ,'> overwritten.']);
            end

            %construct dataset
            topr = [obj.s, obj.q, obj.mn];
            
            %open file with option to append data
            fid = fopen(obj.fname, 'a');
            
            %print each line
            for k = 1:numel(obj.s)
                fprintf(fid, '\n');
                fprintf(fid, obj.FORMAT, topr(k, :));
            end
            
            %close file
            fclose(fid);
        end
        
    end
end