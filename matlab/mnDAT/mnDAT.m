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
        fname = {}; %name of file with path
        
        kilca;      %boolean: KiLCA type or not
        
        n           %toroidal modenumber associated to this file
        
        s           %flux surface label: psi or r for KiLCA
        q           %safety factor
        mn          %mn data
    end
    
    methods
        function obj = mnDAT(path, q, i, n, kilca)
            %##############################################################
            %function obj = mnDAT(path, q, i, n, kilca)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % creates an instance of the class.
            %##############################################################
            % path  ... path of the file
            % q     ... quantitiy contained in mn. e.g.: B
            % i     ... index of q. e.g.: r (-> B_r)
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
            obj.fname = [p, '/', q, 'm', num2str(n), '_', i, k, '.dat'];
            
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
                error('mn must be odd in columns.');
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
                 
            %check if file exists and delete
            if isfile(obj.fname)
                delete(obj.fname);
                warning(['existing file <', obj.fname ,'> overwritten.']);
            end

            %construct dataset
            topr = [obj.s, obj.q, obj.mn];
            
            %open file with option to append data
            fid = fopen(obj.fname, 'a');
            
            %print header line
            %head = {'psi', 'ni[m^-3]', 'ne[m^-3]', 'ti[eV]', 'te[eV]', 'we[rad/s]'};
            %cellfun(@(s) fprintf(fid, '%16s ' , s), head);
            
            %print each line
            for k = 1:numel(obj.s)
                fmt = '%16.10e ';
                fprintf(fid, '\n');
                fprintf(fid, fmt, topr(k, :));
            end
            
            %close file
            fclose(fid);
        end
        
    end
end