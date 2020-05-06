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
    end
    
    properties (Access=private)
        r_priv
        de11_priv
        de12_priv
        de22_priv
        di11_priv
        di12_priv
        di22_priv
        Br_Abs_priv
        Br_minus_ckpEs_over_wexB_priv
        Br_minus_cksEp_over_wexB_priv
        Je_Abs_priv
        Ji_Abs_priv
        Je_plus_Ji_Abs_priv
    end
    
    properties (Dependent)
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
    
    %get for dependent properties
    methods
       function q = get.r(obj)
           if(isempty(obj.r_priv))
               raw = load(obj.path);
               obj.r_priv = raw(:, 1);
           end
           q = obj.r_priv;
       end
       function q = get.de11(obj)
           if(isempty(obj.de11_priv))
               raw = load(obj.path);
               obj.de11_priv = raw(:, 2);
           end
           q = obj.de11_priv;
       end
       function q = get.de12(obj)
           if(isempty(obj.de12_priv))
               raw = load(obj.path);
               obj.de12_priv = raw(:, 3);
           end
           q = obj.de12_priv;
       end
       function q = get.de22(obj)
           if(isempty(obj.de22_priv))
               raw = load(obj.path);
               obj.de22_priv = raw(:, 4);
           end
           q = obj.de22_priv;
       end
       function q = get.di11(obj)
           if(isempty(obj.di11_priv))
               raw = load(obj.path);
               obj.di11_priv = raw(:, 5);
           end
           q = obj.di11_priv;
       end
       function q = get.di12(obj)
           if(isempty(obj.di12_priv))
               raw = load(obj.path);
               obj.di12_priv = raw(:, 6);
           end
           q = obj.di12_priv;
       end
       function q = get.di22(obj)
           if(isempty(obj.di22_priv))
               raw = load(obj.path);
               obj.di22_priv = raw(:, 7);
           end
           q = obj.di22_priv;
       end
       function q = get.Br_Abs(obj)
           if(isempty(obj.Br_Abs_priv))
               raw = load(obj.path);
               obj.Br_Abs_priv = raw(:, 8);
           end
           q = obj.Br_Abs_priv;
       end
       function q = get.Br_minus_ckpEs_over_wexB(obj)
           if(isempty(obj.Br_minus_ckpEs_over_wexB_priv))
               raw = load(obj.path);
               obj.Br_minus_ckpEs_over_wexB_priv = raw(:, 9);
           end
           q = obj.Br_minus_ckpEs_over_wexB_priv;
       end
       function q = get.Br_minus_cksEp_over_wexB(obj)
           if(isempty(obj.Br_minus_cksEp_over_wexB_priv))
               raw = load(obj.path);
               obj.Br_minus_cksEp_over_wexB_priv = raw(:, 10);
           end
           q = obj.Br_minus_cksEp_over_wexB_priv;
       end
       function q = get.Je_Abs(obj)
           if(isempty(obj.Je_Abs_priv))
               raw = load(obj.path);
               obj.Je_Abs_priv = raw(:, 11);
           end
           q = obj.Je_Abs_priv;
       end
       function q = get.Ji_Abs(obj)
           if(isempty(obj.Ji_Abs_priv))
               raw = load(obj.path);
               obj.Ji_Abs_priv = raw(:, 12);
           end
           q = obj.Ji_Abs_priv;
       end
       function q = get.Je_plus_Ji_Abs(obj)
           if(isempty(obj.Je_plus_Ji_Abs_priv))
               raw = load(obj.path);
               obj.Je_plus_Ji_Abs_priv = raw(:, 13);
           end
           q = obj.Je_plus_Ji_Abs_priv;
       end
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
        end
    end
end
