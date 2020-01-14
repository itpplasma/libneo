classdef KiLCA_multimode < handle
%classdef KiLCA_multimode < handle
%##########################################################################
% description of class:
%--------------------------------------------------------------------------
% This class is used to combine the postprocessing of multiple modes and
% make figures for quantities over all modes.
%##########################################################################
% properties:
%--------------------------------------------------------------------------
% *) col, row
% *) postprocessors
%##########################################################################
% methods:
%--------------------------------------------------------------------------
% *) function obj = KiLCA_multimode(post)
% *) function plotResWidth(obj, varargin)
% *) function plotQuant(obj, quant, label, type, varargin)
% *) function plotVsVac(obj, vac, quant, label, type, varargin)
% *) function plotShielding(obj, vac, varargin)
%##########################################################################

    %author:   Philipp Ulbl
    %created:  xx.11.2019
        
    properties (Dependent=true)
        col     %dynamically calculates the optimal number of columns for multimode plots
        row     %dynamically calculates the optimal number of rows for multimode plots
    end
    
    methods
        function q = get.col(obj)
            %optimal number of columns for multimode plots
            q = ceil(sqrt(3 * numel(obj.postprocessors) / 4));
        end 
        function q = get.row(obj)
            %optimal number of rows for multimode plots
            q = ceil(numel(obj.postprocessors) / obj.col);
        end
    end
    
    properties
        postprocessors
    end
    
    methods
        function obj = KiLCA_multimode(post)
            %##############################################################
            %function obj = KiLCA_multimode(post)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % Constructor of the KiLCA_multimode class. Needs all
            % postprocessors from KiLCA_interface parent as a cell array of
            % KiLCA_postprocessor class.
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % pos   ... KiLCA_postprocessors
            %##############################################################
            
            obj.postprocessors = post;
        end
        
        function plotResWidth(obj, varargin)
            %##############################################################
            %function plotResWidth(obj, varargin)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % Creates a plot with a subplot for each mode containing the
            % parallel current, the location of the resonant surface and 2
            % lines indicating the width of the resonant surface.
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % varargin ... plot arguments
            %##############################################################
            
            figure('units', 'normalized', 'outerposition', [0, 0, 1, 1]);
            axis tight

            for k=1:numel(obj.postprocessors)
                post = obj.postprocessors{k};
                subplot(obj.row, obj.col, k)
                %plot J parallel (res surface included)
                post.plot_single('Jpar', 'J_{||} / statA cm^{-2}', 'Abs', varargin{:})
                hold on
                %plot res with lines
                plot([post.rres-post.d, post.rres-post.d], ylim, ':r', 'LineWidth', 2);
                plot([post.rres+post.d, post.rres+post.d], ylim, ':r', 'LineWidth', 2);
                hold off

                title([vec2str(post.mode, '%d'), ', d = ', num2str(post.d)])
            end
        end
        
        function plotQuant(obj, quant, label, type, varargin)
            %##############################################################
            %function plotQuant(obj, quant, label, type, varargin)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % Creates a plot with any quantity plotted for each mode in a
            % subplot. The corresponding location of the resonant surface
            % is included.
            %##############################################################
            % input:
            %--------------------------------------------------------------     
            % quant    ... quantity to plot: string matching a property in 
            %              KiLCA_postprocessor class.
            % label    ... ylabel of the plot       
            % type     ... type: Abs, Re or Im      
            % varargin ... plot arguments
            %##############################################################
            
            figure('units', 'normalized', 'outerposition', [0, 0, 1, 1]);
            axis tight

            for k=1:numel(obj.postprocessors)
                post = obj.postprocessors{k};
                subplot(obj.row, obj.col, k)
                post.plot_single(quant, label, type, varargin{:})
                ax = gca;
                title([ax.Title.String, ' for ', vec2str(post.mode, '%d')])
            end
        end
        
        function plotVsVac(obj, vac, quant, label, type, varargin)
            %##############################################################
            %function plotVsVac(obj, vac, quant, label, type, varargin)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % Creates a plot with any quantity plotted for each mode in a
            % subplot vs. the same quantity for a KiLCA vacuum run.
            % Requires obj to be a flre KiLCA interface.
            % The corresponding location of the resonant surface
            % is included.
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % vac      ... KiLCA_interface for vacuum run
            % quant    ... quantity to plot: string matching a property in 
            %              KiLCA_postprocessor class.
            % label    ... ylabel of the plot       
            % type     ... type: Abs, Re or Im      
            % varargin ... plot arguments
            %##############################################################
            
            %check type of vac and obj
            if(~strcmp(vac.run_type, 'vacuum'))
                error('vac is not a vacuum run.')
            end
            if(~strcmp(obj.run_type, 'flre'))
                error('obj is not a flre run.')
            end
            
            %plot for vac first
            vac.multimode.plotQuant(quant, label, type, 'DisplayName', 'vacuum', varargin{:})

            %include plots for flre
            for k=1:numel(obj.postprocessors)
                post = obj.postprocessors{k};
                subplot(obj.row, obj.col, k)
                hold on
                post.plot_single(quant, label, type, 'DisplayName', 'plasma', varargin{:})
                hold off
                ax = gca;
                title([ax.Title.String, ' for ', vec2str(post.mode, '%d')])
            end
            
            legend('vacuum', 'q=m/n', 'plasma')
            legend('Location', 'NorthWest')
        end
        
        function plotShielding(obj, vac, varargin)
            %##############################################################
            %function plotShielding(obj, vac, varargin)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % Creates a plot with the formfactors at the resonant surface
            % location at their corresponding position on r.
            % Requires obj to be a flre KiLCA interface.
            %##############################################################
            % input:
            %--------------------------------------------------------------
            % vac      ... KiLCA_interface for vacuum run   
            % varargin ... plot arguments
            %##############################################################
            
            %check type of vac and obj
            if(~strcmp(vac.run_type, 'vacuum'))
                error('vac is not a vacuum run.')
            end
            if(~strcmp(obj.run_type, 'flre'))
                error('obj is not a flre run.')
            end
            
            figure('units', 'normalized', 'outerposition', [0, 0, 1, 0.65]);
            axis tight
            
            %initialize res pos and form factor (ratio)
            res = zeros(size(obj.postprocessors));
            ratio = zeros(size(obj.postprocessors));

            hold on
            %read res and calc ratio for each mode
            for k=1:numel(obj.postprocessors)
                
                %extract postprocessors
                post = obj.postprocessors{k};
                postv = vac.multimode.postprocessors{k};
                
                res(k) = post.rres;
                ratio(k) = interp1(post.r, abs(post.Br), post.rres) / interp1(postv.r, abs(postv.Br), post.rres);
                
                %include text with m,n
                text(res(k) + sum(xlim .* [0, 0.01]), ratio(k) + sum(ylim .* [0, 0.001]), vec2str(post.mode, '%d'))
            end
            
            plot(res, ratio, '-ob', 'MarkerSize', 5, 'MarkerFaceColor', 'b');
            hold off

            xlim([post.rmin, post.ra])    
            ylim([0, min(1, 2 * max(ratio))])
            
            xlabel('r / cm')
            ylabel('|B^{plas}_{r}| / |B^{vac}_{r}| at q = m/n')
            
            %title('Form factors for different modenumbers at q = res')
        end
    end
end