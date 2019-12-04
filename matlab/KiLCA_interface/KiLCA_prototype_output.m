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
%modified: 04.11.2019
    
    properties
    end
    
    methods (Access = 'public')
        
    end
    
    methods (Access = 'public')
        function p = plot_single(obj, a, u, type, varargin)
            %##############################################################
            %function p = plot_single(obj, a, u, type, varargin)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % plots single property over radius given by name a
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % a         ... property to be plot
            % u         ... ylabel as text
            % type      ... type of plot: Re, Im or Abs (=empty)
            % varargin  ... plot arguments (if used, type must be non-empty!)
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % p         ... plot handle
            %##############################################################
                        
            dim = size(obj.(a));
            if(dim(2) > 1)
                x = obj.(a)(:, 1);
                y = obj.(a)(:, 2);
            else
                x = obj.r;
                y = obj.(a);
            end
            
            if(strcmp(type, 'Re'))
                y = real(y);
            elseif(strcmp(type, 'Im'))
                y = imag(y);
            elseif(strcmp(type, 'Abs'))
                y = abs(y);
            end
            
            p = plot(x, y, 'DisplayName', u, varargin{:});
            title(obj.format_title(a, type))
            xlabel('r / cm')
            ylabel(u)
        end
        
        function plot_triple(obj, a, b, c, u, type, varargin)
            %##############################################################
            %function p = plot_single(obj, a, u, type, varargin)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % plots a set of 3 properties over radius given by a,b,c
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % a, b, c   ... properties to be plot
            % u         ... ylabel as text
            % type      ... type of plot: Re, Im or Abs (=empty)
            % varargin  ... plot arguments (if used, type must be non-empty!)
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % p         ... plot handle
            %##############################################################
            
            subplot(1, 3, 1)
            plot_single(obj, a, u, type, varargin{:});
            
            subplot(1, 3, 2)
            plot_single(obj, b, u, type, varargin{:});
            
            subplot(1, 3, 3)
            plot_single(obj, c, u, type, varargin{:});
        end
    end
    
    methods (Static, Access = 'private') 
        function a = format_title(a, type)
            %converts property names to appropriate plot titles
            
            
            if contains(a, '0')
                a = strrep(a, '0', '_0');
                if contains(a, 'b')
                    a = strrep(a, 'b', 'B');
                end
                if contains(a, 'j')
                    a = strrep(a, 'j', 'J');
                end
            end
            
            if contains(a, 'par')
                a = strrep(a, 'par', '_{\mid\mid}');
            elseif contains(a, 'perp')
                a = strrep(a, 'perp', '_{\perp}');
            %r underscore
            elseif contains(a, 'r')
                a = strrep(a, 'r', '_r');
            %theta underscore
            elseif contains(a, 'th')
                a = strrep(a, 'th', '_{\theta}');
            %z underscore
            elseif contains(a, 'z')
                a = strrep(a, 'z', '_z');
            end
            
            %absolute symbol
            if contains(a, '_Abs') || strcmp(type, 'Abs')
                a = erase(a, '_Abs');
                if(~contains(a, '0'))
                    a = ['|', a, '|'];
                end
            %real part
            elseif contains(a, '_Re') || strcmp(type, 'Re')
                a = erase(a, '_Re');
                if(~contains(a, '0'))
                    a = ['Re(', a, ')'];
                end
            %imag part
            elseif contains(a, '_Im') || strcmp(type, 'Im')
                a = erase(a, '_Im');
                if(~contains(a, '0'))
                    a = ['Im(', a, ')'];
                end
            end
            
            %erase not needed index
            if contains(a, '_m')
                a = erase(a, '_m');
            elseif contains(a, '_p')
                a = erase(a, '_p');
            end
            
        end
    end
end