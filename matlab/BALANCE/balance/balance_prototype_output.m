classdef (Abstract) balance_prototype_output < handle
%##########################################################################
%classdef (Abstract) balance_prototype_output < handle
%##########################################################################
% description of class:
%--------------------------------------------------------------------------
% This class is a prototype for output classes of the Balance code.
%##########################################################################
% properties:
%--------------------------------------------------------------------------
%##########################################################################
% methods:
%--------------------------------------------------------------------------
% *) p = plot(obj, a, u, varargin)
%########################################################################## 

%author:   Philipp Ulbl
%created:  13.01.2020

    properties
        
    end
    
    methods
        function p = plot(obj, a, u, varargin)
            %##############################################################
            %function p = plot(obj, a, u, type, varargin)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % plots property over radius given by name a
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % a         ... property to be plot
            % u         ... ylabel as text
            % varargin  ... plot arguments (if used, type must be non-empty!)
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % p         ... plot handle
            %##############################################################
                        
            x = obj.r;
            y = obj.(a);
            
            p = plot(x, y, 'DisplayName', u, varargin{:});
            title(obj.format_title(a))
            xlabel('r / cm')
            ylabel(u)
        end
    end
        
    methods (Static, Access = 'private') 
        function a = format_title(a)
            %converts property names to appropriate plot titles
            
            %absolute symbol
            if contains(a, '_Abs')
                a = erase(a, '_Abs');
                a = ['|', a, '|'];
            end
            
            %remove underlines
            if contains(a, '_')
                a = strrep(a, '_', ' ');
            end
            
            %minus sign
            if contains(a, 'minus')
                a = strrep(a, 'minus', '-');
            end
            
            %plus sign
            if contains(a, 'plus')
                a = strrep(a, 'plus', '+');
            end
            
            %over /
            if contains(a, 'over')
                a = strrep(a, 'over', '/');
            end
            
            %square root
            if contains(a, 'Sqrt')
                a = strrep(a, 'Sqrt', '\surd');
            end
            
            %theta
            if contains(a, 'theta')
                a = strrep(a, 'theta', '_\theta');
            end
        end
    end
end