classdef f1000 < balance_prototype_output
%##########################################################################
%classdef f1000 < output
%##########################################################################
% description of class:
%--------------------------------------------------------------------------
% This class is used to read fort.1000 files which are output of the
% Balance code.
%##########################################################################
% properties:
%--------------------------------------------------------------------------
% *) path
% *) r, n, Vz, Te, Ti, Er, Sqrtg_Btheta_over_c
%##########################################################################
% methods:
%--------------------------------------------------------------------------
% *) obj = f1000(fpath)
%########################################################################## 

%author:   Philipp Ulbl
%created:  13.01.2020

    properties
        path
        
        r                   %effective radius
        n                   %density
        Vz                  %toroidal velocity
        Te                  %electron temperature
        Ti                  %ion temperature
        Er                  %radial electric field
        Sqrtg_Btheta_over_c %
    end
    
    methods
        function obj = f1000(fpath)
            %##############################################################
            %function obj = f1000(fpath)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % constructor of the f1000 class.
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % path  ... path of input file
            %############################################################## 
            
            obj.path = fpath;
            
            %load data
            raw = load(fpath);
            
            %assign columns to properties
            obj.r = raw(:, 1);
            obj.n = raw(:, 2);
            obj.Vz = raw(:, 3);
            obj.Te = raw(:, 4);
            obj.Ti = raw(:, 5);
            obj.Er = raw(:, 6);
            obj.Sqrtg_Btheta_over_c = raw(:, 7);
        end
    end
end
