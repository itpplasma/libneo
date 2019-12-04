classdef KiLCA_multimode
    %KILCA_MULTIMODE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Dependent=true)
        col
        row
    end
    
    methods
        function q = get.col(obj)
           q = ceil(sqrt(3 * numel(obj.postprocessors) / 4));
        end 
        function q = get.row(obj)
           q = ceil(numel(obj.postprocessors) / obj.col);
        end
    end
    
    properties
        postprocessors
    end
    
    methods
        function obj = KiLCA_multimode(post)
            
            obj.postprocessors = post;
        end
        
        function plotResWidth(obj, varargin)
            
            figure('units', 'normalized', 'outerposition', [0, 0, 1, 1]);
            axis tight

            for k=1:numel(obj.postprocessors)
                post = obj.postprocessors{k};
                subplot(obj.row, obj.col, k)
                post.plot_single('Jpar', 'J_{||} / statA cm^{-2}', 'Abs', varargin{:})
                hold on
                plot([post.rres-post.d, post.rres-post.d], ylim, ':r', 'LineWidth', 2);
                plot([post.rres+post.d, post.rres+post.d], ylim, ':r', 'LineWidth', 2);
                hold off

                title([vec2str(post.mode, '%d'), ', d = ', num2str(post.d)])
            end
            
        end
        
        function plotQuant(obj, quant, label, type, varargin)
            
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
            
            vac.multimode.plotQuant(quant, label, type, 'DisplayName', 'vacuum', varargin{:})

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
                
            figure('units', 'normalized', 'outerposition', [0, 0, 1, 0.65]);
            axis tight
            
            res = zeros(size(obj.postprocessors));
            ratio = zeros(size(obj.postprocessors));

            hold on
            for k=1:numel(obj.postprocessors)
                
                post = obj.postprocessors{k};
                postv = vac.multimode.postprocessors{k};
                
                res(k) = post.rres;
                ratio(k) = interp1(post.r, abs(post.Br), post.rres) / interp1(postv.r, abs(postv.Br), post.rres);
    
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

