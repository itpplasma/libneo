classdef KiLCA_modes < KiLCA_prototype_input
%classdef KiLCA_modes
%##########################################################################
% description of class:
%--------------------------------------------------------------------------
% class containing the information of the antenna modes in KiLCA.
%##########################################################################
% properties:
%--------------------------------------------------------------------------
% *) m, n; default = 3, 2
% READONLY:
% *) INDICES, BLUEPRINT, READY
%##########################################################################
% methods:
%--------------------------------------------------------------------------
% *) function obj = KiLCA_modes(m, n)
% *) function c = plain(obj)
%##########################################################################

%author:   Philipp Ulbl
%created:  08.08.2019
%modified: 21.08.2019

    properties (SetAccess = 'private')
        BLUEPRINT = 'modes.in';   %name of blueprint file
        READY = true;             %flag: ready to run
        
        m = 3
        n = 2
    end
    
    methods
        function obj = KiLCA_modes(m, n)
            %##############################################################
            %function obj = KiLCA_modes(m, n)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % constructor of the modes class. optional one can set m,n
            % initially.
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % m     ... poloidal modenumber (optional, default = 3)
            % n     ... toroidal modenumber (optional, default = 2)
            %##############################################################
            
            if nargin == 2
                obj.set(m, n);
            elseif nargin == 1
                obj.set(m, []);
            end
        end
        
        function set(obj, m, n)
            %##############################################################
            %function set(obj, m, n)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % set one of the modenumbers or both
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % m     ... poloidal modenumber (optional, default = 3)
            % n     ... toroidal modenumber (optional, default = 2)
            %##############################################################
            
            if ~isempty(m)
                if all(size(m) > 1)
                    error('m must be a scalar or vector.')
                end
                obj.m = m;
            end
            
            if nargin > 2 && ~isempty(n)
                if all(size(n) > 1)
                    error('n must be a scalar or vector.')
                end
                if ~all(size(m) == size(n))
                    error('size of m and n should be the same.')
                end
                obj.n = n;
            end
        end
        
        function c = plain(obj)
            %##############################################################
            %function c = plain(obj)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % returns the class as a cell-array with the name of the
            % properties in each cell. properties which are structs will be
            % split into their elements.
            % *overrides superclass .plain() method
            %##############################################################
            % output:
            %--------------------------------------------------------------
            % c     ... cell array
            %##############################################################
            
            c = cell(1, numel(obj.m));
            
            for k = 1:numel(obj.m)
                c{k} = vec2str([obj.m(k), obj.n(k)], '%d');
            end
        end
        function write(obj, ~, path_to)
            %##############################################################
            %function write(obj, path_from, path_to)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % writes properties of the class into a KiLCA input file.
            % *overrides superclass .write() method
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % path_from  ... path where the blueprint is from
            % path_to    ... path where the input file will be written
            %##############################################################
            
            save_file(obj.plain(), [path_to, obj.BLUEPRINT]);
        end
    end
end