classdef KiLCA_modes < handle & blueprint
%classdef KiLCA_modes < handle & blueprint
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
% *) function gen_modes(obj, mlim, nlim, qlim)
% *) function write(obj, ~, path_to)
%##########################################################################

%author:   Philipp Ulbl
%created:  08.08.2019
%modified: 26.03.2020

    properties (Transient, SetAccess = 'protected')
        INDICES = [];
        BLUEPRINT = 'modes.in';   %name of blueprint file
        SEP = '#'
    end
    properties (SetAccess = 'private')
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
            obj.READY = true;
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
        
        function gen_modes(obj, mlim, nlim, qlim)
            %##############################################################
            %function gen_modes(obj, mlim, nlim, qlim)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % generates all possible unique mode combinations of m between
            % mlim(1) and mlim(2) and n between nlim(1) and nlim(2) with
            % m/n greater than 1 and smaller than qlim. Sorts with
            % ascending m/n. Sets m and n property of the class.
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % mlim  ... 2-entry vector with min(m) and max(m)
            % nlim  ... 2-entry vector with min(n) and max(n)
            % qlim  ... maximum q (abs)
            %##############################################################
                       
            %initialize vectors from lim(1) to lim(2)
            m1 = mlim(1):mlim(2);
            n1 = nlim(1):nlim(2);
            
            %generate all possible combinations and map to vectors again
            [M, N] = meshgrid(m1, n1);
            m1 = M(:);
            n1 = N(:);
            
            %filter out all m <= n (because q always > 1)
            I = find(m1<=n1);
            m1(I) = [];
            n1(I) = [];

            %filter out repeating m/n combinations:
            %e.g. (4,2) and (8,4)
            [~, I] = unique(m1./n1);
            m1 = m1(I);
            n1 = n1(I);

            %filter out all m/n combinations that are larger than qlim
            I = m1./n1 <= qlim;
            m1 = m1(I);
            n1 = n1(I);

            %sort according to m/n ascending
            [~, I] = sort(m1./n1);
            m1 = m1(I);
            n1 = n1(I);

            %set modes
            obj.set(m1, n1)
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