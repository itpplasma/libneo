classdef (Abstract) KiLCA_prototype_output < handle
%classdef (Abstract) KiLCA_prototype_output < handle
%##########################################################################
% description of class:
%--------------------------------------------------------------------------
% prototype class for all KiLCA output classes.
%##########################################################################
% properties:
%--------------------------------------------------------------------------
%##########################################################################
% methods:
%--------------------------------------------------------------------------
% *) 
%##########################################################################

%author:   Philipp Ulbl
%created:  21.08.2019
%modified: 21.08.2019
    
    properties
    end
    
    methods (Access = 'public')
        
    end
    
    methods (Access = 'public')
        
        function p = plot_single(obj, a, u, varargin)
            %plots single property over radius given by name a
            
            p = plot(obj.(a)(:, 1), obj.(a)(:, 2), varargin{:});
            title(obj.format_title(a))
            xlabel('r / cm')
            ylabel(u)
        end
        
        function plot_triple(obj, a, b, c, u, varargin)
            %plots a set of 3 properties over radius given by a,b,c
            
            subplot(1, 3, 1)
            plot_single(obj, a, u, varargin{:});
            
            subplot(1, 3, 2)
            plot_single(obj, b, u, varargin{:});
            
            subplot(1, 3, 3)
            plot_single(obj, c, u, varargin{:});
        end
    end
    
    methods (Static, Access = 'private') 
        function a = format_title(a)
            %converts property names to appropriate plot titles
            
            %r underscore
            if contains(a, 'r')
                a = strrep(a, 'r', '_r');
            %theta underscore
            elseif contains(a, 'th')
                a = strrep(a, 'th', '_{\theta}');
            %z underscore
            elseif contains(a, 'z')
                a = strrep(a, 'z', '_z');
            end
            
            %absolute symbol
            if contains(a, '_Abs')
                s = erase(a, '_Abs');
                a = ['|', s, '|'];
            %real part
            elseif contains(a, '_Re')
                s = erase(a, '_Re');
                a = ['Re(', s, ')'];
            %imag part
            elseif contains(a, '_Im')
                s = erase(a, '_Im');
                a = ['Im(', s, ')'];
            end
        end
    end
end