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
           q = obj.saveGetProp('r_priv');
       end
       function q = get.de11(obj)
           q = obj.saveGetProp('de11_priv');
       end
       function q = get.de12(obj)
           q = obj.saveGetProp('de12_priv');
       end
       function q = get.de22(obj)
           q = obj.saveGetProp('de22_priv');
       end
       function q = get.di11(obj)
           q = obj.saveGetProp('di11_priv');
       end
       function q = get.di12(obj)
           q = obj.saveGetProp('di12_priv');
       end
       function q = get.di22(obj)
           q = obj.saveGetProp('di22_priv');
       end
       function q = get.Br_Abs(obj)
           q = obj.saveGetProp('Br_Abs_priv');
       end
       function q = get.Br_minus_ckpEs_over_wexB(obj)
           q = obj.saveGetProp('Br_minus_ckpEs_over_wexB_priv');
       end
       function q = get.Br_minus_cksEp_over_wexB(obj)
           q = obj.saveGetProp('Br_minus_cksEp_over_wexB_priv');
       end
       function q = get.Je_Abs(obj)
           q = obj.saveGetProp('Je_Abs_priv');
       end
       function q = get.Ji_Abs(obj)
           q = obj.saveGetProp('Ji_Abs_priv');
       end
       function q = get.Je_plus_Ji_Abs(obj)
           q = obj.saveGetProp('Je_plus_Ji_Abs_priv');
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
        
        function b = getLastBrAbs(obj)
            %##############################################################
            %function b = getLastBrAbs(obj)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % returns the last entry of Br_Abs in an efficient way.
            %##############################################################
            % output:
            %--------------------------------------------------------------
            % b     ... last entry of Br_Abs
            %############################################################## 
            
            %get last line using system command
            [~, raw] = system(['tail -n 1 ', obj.path]);
            %split by space deliminator
            raw = strsplit(raw, ' ');
            %remove first and last entry (no chars)
            raw = raw(2:(end-1));
            %return value in 8th entry(Br_Abs)
            b=str2double(raw{8});
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
            raw = reshape(raw{1}, 13, numel(raw{:})/13)';
            %assign columns to properties
            obj.r_priv = raw(:, 1);
            obj.de11_priv = raw(:, 2);
            obj.de12_priv = raw(:, 3);
            obj.de22_priv = raw(:, 4);
            obj.di11_priv = raw(:, 5);
            obj.di12_priv = raw(:, 6);
            obj.di22_priv = raw(:, 7);
            obj.Br_Abs_priv = raw(:, 8);
            obj.Br_minus_ckpEs_over_wexB_priv = raw(:, 9);
            obj.Br_minus_cksEp_over_wexB_priv = raw(:, 10);
            obj.Je_Abs_priv = raw(:, 11);
            obj.Ji_Abs_priv = raw(:, 12);
            obj.Je_plus_Ji_Abs_priv = raw(:, 13);
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
