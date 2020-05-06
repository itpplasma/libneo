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
    end
    
    properties (Access=private)
        r_priv
        n_priv
        Vz_priv
        Te_priv
        Ti_priv
        Er_priv
        Sqrtg_Btheta_over_c_priv
    end
    
    properties (Dependent)
        r                   %effective radius
        n                   %density
        Vz                  %toroidal velocity
        Te                  %electron temperature
        Ti                  %ion temperature
        Er                  %radial electric field
        Sqrtg_Btheta_over_c %
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
       function q = get.n(obj)
           if(isempty(obj.n_priv))
               raw = load(obj.path);
               obj.n_priv = raw(:, 2);
           end
           q = obj.n_priv;
       end
       function q = get.Vz(obj)
           if(isempty(obj.Vz_priv))
               raw = load(obj.path);
               obj.Vz_priv = raw(:, 3);
           end
           q = obj.Vz_priv;
       end
       function q = get.Te(obj)
           if(isempty(obj.Te_priv))
               raw = load(obj.path);
               obj.Te_priv = raw(:, 4);
           end
           q = obj.Te_priv;
       end
       function q = get.Ti(obj)
           if(isempty(obj.Ti_priv))
               raw = load(obj.path);
               obj.Ti_priv = raw(:, 5);
           end
           q = obj.Ti_priv;
       end
       function q = get.Er(obj)
           if(isempty(obj.Er_priv))
               raw = load(obj.path);
               obj.Er_priv = raw(:, 6);
           end
           q = obj.Er_priv;
       end
       function q = get.Sqrtg_Btheta_over_c(obj)
           if(isempty(obj.Sqrtg_Btheta_over_c_priv))
               raw = load(obj.path);
               obj.Sqrtg_Btheta_over_c_priv = raw(:, 7);
           end
           q = obj.Sqrtg_Btheta_over_c_priv;
       end
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
        end
    end
end