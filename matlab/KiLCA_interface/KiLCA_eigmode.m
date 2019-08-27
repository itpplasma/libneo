classdef KiLCA_eigmode < KiLCA_prototype_input
%classdef KiLCA_eigmode
%##########################################################################
% description of class:
%--------------------------------------------------------------------------
% class containing the information of the eigmode in KiLCA.
%##########################################################################
% properties:
%--------------------------------------------------------------------------
% *) output, flag_fscan
% *) fgrid_redim, fgrid_remin, fgrid_remax, fgrid_imdim, fgrid_immin,
%    fgrid_immax
% *) flag_stopcrit, det_abs, rootseq_abs, rootseq_rel 
% *) df
% *) flag_testroots, flag_deb
% *) rsearch_nstart, rsearch_istart, rsearch_iend, rsearch_points
% READONLY:
% *) INDICES, BLUEPRINT, READY
%##########################################################################
% methods:
%--------------------------------------------------------------------------
% *) function obj = KiLCA_eigmode()
% *) function c = plain(obj)
%##########################################################################

%author:   Philipp Ulbl
%created:  21.08.2019
%modified: 21.08.2019

    properties (SetAccess = 'private')
        INDICES = [2, 5, 8:13, 16:19, 22, 25:26, 29:31, 34:40];  %indices of parameters in blueprint files
        BLUEPRINT = 'eigmode.in';                                %name of blueprint file
        READY = true;                                            %flag: ready to run
    end
    
    properties
        output = 'roots.dat'    %file name with determinat values (either from roots search, or on a grid)
        flag_fscan = 1          %choose if to use frequency scan or roots search: 0 - if roots search is used, 1 - if det is evaluated on a frequency grid
        
        %frequency grid settings:
        
        fgrid_redim = 1         %dimension of a grid for real part
        fgrid_remin = 0         %minimum real value of frequency
        fgrid_remax = 1e6       %maximim real value of frequency
        fgrid_imdim = 100       %dimension of a grid for imag part
        fgrid_immin = 0         %minimum imag value of frequency
        fgrid_immax = 2.1e5     %maximim imag value of frequency
        
        %Stopping criteria for roots search:
        
        flag_stopcrit = 1       %stopping criteria: 0 - det residial, 1 - root convergence
        det_abs = 1e-14         %absolute value of the determinant
        rootseq_abserr = 1e-10  %absolute error for roots sequence
        rootseq_relerr = 1e-10  %relative error for roots sequence
        
        %For omega derivative:
        
        df = 1e-3               %delta freq used to compute derivative over omega numerically
        
        %Test roots:
        
        flag_testroots = 0      %1 - test roots by contour integration, 0 - skip test
        flag_deb = 0            %flag for debugging
        
        %Starting points indices for roots search:
        rsearch_nstart = 4      %number of starting points
        rsearch_istart = 0      %starting index (from 0)
        rsearch_iend = 3        %ending index   (from 0)
        rsearch_points = [zeros(1, 7); logspace(2, 8, 7)]' %Starting points array: default: [1e1, 1e2, ... 1e7]
    end
    
    methods
        function obj = KiLCA_eigmode()
            %##############################################################
            %function obj = KiLCA_eigmode()
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % constructor of the eigmode class. empty.
            %##############################################################
        end
        
        function c = plain(obj)
            %##############################################################
            %function c = plain(obj)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % returns the class as an cell-array with the name of the
            % properties in each cell. properties which are structs will be
            % split into their elements.
            % *extends superclass .plain() method
            %##############################################################
            % output:
            %--------------------------------------------------------------
            % c     ... cell array
            %##############################################################
            
            %use superclass method
            c = plain@KiLCA_prototype_input(obj);
            
            %modify last cell: split array into cell for each row
            
            %delete last row
            c(end) = [];
            %add other rows
            for k=1:size(obj.rsearch_points, 1)
                 c{end+1} = obj.rsearch_points(k, :);
            end
        end
        
        function write(obj, path_from, path_to)
            %##############################################################
            %function write(obj, path_from, path_to)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % writes properties of the class into a KiLCA input file.
            % *extends superclass .write() method
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % path_from  ... path where the blueprint is from
            % path_to    ... path where the input file will be written
            %##############################################################
            
            %if more than 7 rsearchpoints are used, adjust the INDICES to
            %write them into the input file
            if numel(obj.rsearch_points) > 7
                obj.INDICES = [obj.INDICES,1:(obj.rsearch_points-7)];
            end
            
            %continue with ordinary write of superclass
            write@KiLCA_prototype_input(obj, path_from, path_to);
        end
    end
end