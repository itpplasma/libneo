% Extract data from an axes and save it to file.
%
% input:
% ------
% ax: axes or line, from which to extract x and y data.
% filename: string, name under which to save output M
% nmax: optional integer, number of data elements (lines)? not used if
%   ax is a line [number of children if ax is axes]
% sparsing: optional boolean [false]
% xlimits: optional two element row vector [[-inf, inf]]
%
% output:
% -------
% M:
%
% side effects:
% -------------
% - saves M to file filename in ascii format.
function M = export_pgf(ax, filename, nmax, sparsing, xlimits)

  if nargin < 4 || isempty(sparsing), sparsing=false; end
  if nargin < 5 || isempty(xlimits),  xlimits=[-inf, inf]; end

  npoints_max = 2000;

  if (isa(ax, 'matlab.graphics.axis.Axes'))

    n = numel(ax.Children);
    if nargin < 3 || isempty(nmax),     nmax=n;         end

    xdat = cell(1, nmax);
    ydat = cell(1, nmax);
    for k = 1:nmax
      xdat{k} = ax.Children(k).XData;
      ydat{k} = ax.Children(k).YData;
    end
  else
   % Hopefully this is now a Line 'matlab.graphics.chart.primitive.Line'
   nmax = 1;
   xdat{1} = ax.XData;
   ydat{1} = ax.YData;
  end

  % or arrayfun(@(x) x.XData,ax.Children,'Uni',0)

  npoints = numel(xdat{1});

  M = zeros(npoints, 2*nmax);
  M(:, end-1:-2:1) = reshape(cell2mat(xdat), npoints, nmax);
  M(:, end:-2:2)   = reshape(cell2mat(ydat), npoints, nmax);

  M = M((M(:,1) >= xlimits(1)) & (M(:,1) <= xlimits(2)) | isnan((M(:,1))), :);
  npoints = size(M,1);
  if (sparsing)
    %Goal = 1000 data points
    if (npoints > npoints_max)
      nth = fix(npoints / npoints_max);
      M = M([1:nth:end, end],:);
    end
  end

  save(filename, 'M', '-ascii');
end
