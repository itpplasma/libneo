classdef f5000 < balance_prototype_output
%##########################################################################
%classdef f5000 < output
%##########################################################################
% description of class:
%--------------------------------------------------------------------------
% This class is used to read fort.5000 files which are output of the
% Balance code.
%##########################################################################
% properties:
%--------------------------------------------------------------------------
% *) path
% *) r, de11, de12, de22, di11, di12, di22
% *) Br_Abs, Br_minus_ckpEs_over_wexB, Br_minus_cksEp_over_wexB
% *) Je_Abs, Ji_Abs, Je_plus_Ji_Abs
%##########################################################################
% methods:
%--------------------------------------------------------------------------
% *) obj = f5000(fpath)
%########################################################################## 

%author:   Philipp Ulbl
%created:  13.01.2020

    properties
        path
        
        r       %effective radius
        de11    %transport coefficient for electrons 11
        de12    %transport coefficient for electrons 12
        de22    %transport coefficient for electrons 22
        di11    %transport coefficient for ions 11
        di12    %transport coefficient for ions 12
        di22    %transport coefficient for ions 22
        Br_Abs  %absolute radial magnetic field perturbation
        Br_minus_ckpEs_over_wexB %
        Br_minus_cksEp_over_wexB %
        Je_Abs          %absolute electron current density
        Ji_Abs          %absolute ion current density
        Je_plus_Ji_Abs  %absolute of total current density
    end
    
    methods
        function obj = f5000(fpath)
            %##############################################################
            %function obj = f5000(fpath)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % constructor of the f5000 class.
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
            obj.de11 = raw(:, 2);
            obj.de12 = raw(:, 3);
            obj.de22 = raw(:, 4);
            obj.di11 = raw(:, 5);
            obj.di12 = raw(:, 6);
            obj.di22 = raw(:, 7);
            obj.Br_Abs = raw(:, 8);
            obj.Br_minus_ckpEs_over_wexB = raw(:, 9);
            obj.Br_minus_cksEp_over_wexB = raw(:, 10);
            obj.Je_Abs = raw(:, 11);
            obj.Ji_Abs = raw(:, 12);
            obj.Je_plus_Ji_Abs = raw(:, 13);
        end
    end
end
