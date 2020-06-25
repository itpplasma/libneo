classdef (Sealed, Abstract) ColorManager
%classdef (Sealed, Abstract) ColorManager
%##########################################################################
% description of class:
%--------------------------------------------------------------------------
% This class acts as container for sets of colors for plotting. It can be
% accessed as a static class with a static method that returns a
% colormatrix to be used in plots. Colors sets can be set by their name
% given by GetSetNames().
% 
% Colors are taken from colorbrewer which has nice sets of colors:
% https://colorbrewer2.org
%##########################################################################
% methods (static):
%--------------------------------------------------------------------------
% *) function nam = GetSetNames()
% *) function col = GetColor(num, name)
%##########################################################################

%author:   Philipp Ulbl
%created:  05.06.2020
    
    properties(Constant)
        
        Accent = [127,201,127;
                  190,174,212;
                  253,192,134;
                  255,255,153;
                   56,108,176;
                  240,  2,127;
                  191, 91, 23;
                  102,102,102]./255;
        
        Dark2 = [ 27,158,119;
                 217, 95,  2;
                 117,112,179;
                 231, 41,138;
                 102,166, 30;
                 230,171,  2;
                 166,118, 29;
                 102,102,102]./255;
              
        Set1= [228, 26, 28;
                55,126,184;
                77,175, 74;
               152, 78,163;
               255,127,  0;
               255,255, 51;
               166, 86, 40;
               247,129,191;
               153,153,153]./255;
           
        Set2= [102,194,165;
               252,141, 98;
               141,160,203;
               231,138,195;
               166,216, 84;
               255,217, 47;
               229,196,148;
               179,179,179]./255;
            
        Set3= [141,211,199;
               255,255,179;
               190,186,218;
               251,128,114;
               128,177,211;
               253,180, 98;
               179,222,105;
               252,205,229;
               217,217,217;
               188,128,189;
               204,235,197;
               255,237,111]./255;
            
        Paired =  [166,206,227;
                    31,120,180;
                   178,223,138;
                    51,160, 44;
                   251,154,153;
                   227, 26, 28;
                   253,191,111;
                   255,127,  0;
                   202,178,214;
                   106, 61,154;
                   255,255,153;
                   177, 89, 40]./255;

        Pastel1 = [251,180,174;
                   179,205,227;
                   204,235,197;
                   222,203,228;
                   254,217,166;
                   255,255,204;
                   229,216,189;
                   253,218,236;
                   242,242,242]./255;
               
        Pastel2 = [179,226,205;
                   253,205,172;
                   203,213,232;
                   244,202,228;
                   230,245,201;
                   255,242,174;
                   241,226,204;
                   204,204,204]./255;
    end
    
    methods(Static)
        function nam = GetSetNames()
        %##################################################################
        %function nam = GetSetNames()
        %##################################################################
        % description:
        %------------------------------------------------------------------
        % This functions returns the names of all color sets included in
        % the class.
        %##################################################################
        % output:
        %------------------------------------------------------------------
        % nam       ... cell array with set names
        %##################################################################
        
            nam = properties('ColorManager');
        end
        
        function col = GetColor(num, wrap, name)
        %##################################################################
        %function col = GetColor(num, wrap, name)
        %##################################################################
        % description:
        %------------------------------------------------------------------
        % This functions returns the specified number of colors as a color
        % matrix. The name of the color set can be specified.
        % If num > colors in the set, a warning together with col=[] is
        % returned.
        %##################################################################
        % input:
        %------------------------------------------------------------------
        % num       ... number of color objects to return
        % wrap      ... boolean that indicates to wrap colors if num >
        %               elements in set (optional, default = false)
        % name      ... name of the color set (optional, default = Set1)
        %##################################################################
        % output:
        %------------------------------------------------------------------
        % col       ... color matrix: each row corresponds to the rgb
        %               values for a single color.
        %##################################################################

            %default for name
            if(nargin < 2 || isempty(wrap))
               wrap = false;
            end
            %default for name
            if(nargin < 3 || isempty(name))
               name = 'Set1';
            end
            
            %get colors from set by name
            col = ColorManager.(name);
            
            % if wrap is activated wrap colors around the 1st dimension
            % with use of modulo
            if(~wrap) 
                %get available sets and check if set specified by name is there
                sets = ColorManager.get_available_sets(num);
                ok = any(cellfun(@(x) strcmp(x, name), sets));

                %return colors if ok
                if(ok)
                    ind = 1:num;
                    col = col(ind, :);
                else
                    warning('specified set has not enough colors, empty field returned.')
                    col = [];
                end
            else
                ind = mod((1:num)-1, size(col, 1))+1;
                col = col(ind, :);
            end
               
        end
        
    end
    
    methods (Static, Access=private)
        function av = get_available_sets(num)
            %gets all sets that contain more or equal elements than num
            
            %get all properties (=sets) of class
            prop = properties('ColorManager');
            av = {};
            
            %test if set has more colors than requested, if yes, add to av
            for k=1:numel(prop)
                if(size(ColorManager.(prop{k}), 1) >= num)
                    av = [av; prop{k}];
                end
            end
        end
    end
end