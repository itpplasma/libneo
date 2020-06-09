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
           q = obj.saveGetProp('r_priv');
       end
       function q = get.n(obj)
           q = obj.saveGetProp('n_priv');
       end
       function q = get.Vz(obj)
           q = obj.saveGetProp('Vz_priv');
       end
       function q = get.Te(obj)
           q = obj.saveGetProp('Te_priv');
       end
       function q = get.Ti(obj)
           q = obj.saveGetProp('Ti_priv');
       end
       function q = get.Er(obj)
           q = obj.saveGetProp('Er_priv');
       end
       function q = get.Sqrtg_Btheta_over_c(obj)
           q = obj.saveGetProp('Sqrtg_Btheta_over_c_priv');
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
    
    methods(Access = private)
        
        function loadFile(obj)
            %loads full file into properties
            
            %open file
            fid = fopen(obj.path);
            %textscan with right format
            raw = textscan(fid, '%25n', 'Delimiter' , ' ', 'MultipleDelimsAsOne', true);
            %reshape single column into 7 columns
            raw = reshape(raw{1}, 7, numel(raw{:})/7)';
            %assign columns to properties
            obj.r_priv = raw(:, 1);
            obj.n_priv = raw(:, 2);
            obj.Vz_priv = raw(:, 3);
            obj.Te_priv = raw(:, 4);
            obj.Ti_priv = raw(:, 5);
            obj.Er_priv = raw(:, 6);
            obj.Sqrtg_Btheta_over_c_priv = raw(:, 7);
            %close file
            fclose(fid);
        end    
        function q = saveGetProp(obj, prop)
            %checks if property has been loaded and loads it if not
            
            if(isempty(obj.(prop)))
                obj.loadFile();
            end
            
            q = obj.(prop);
        end
    end
end