%##########################################################################
% script_test_smooth2level.m
%##########################################################################
% description:
%--------------------------------------------------------------------------
% script to test function smooth2level
%##########################################################################
     
%author:   Philipp Ulbl
%created:  21.01.2020

%NOTICE: smoothing to arbitrary level is not a good idea.
%        you can try with Er.dat using level=1 and level=2 or with
%        Te.dat using level=3 and level=4. datapoints themselves may get
%        off the original ones if you use too high level

%smoothing level
level = 2;

%test data
prof = './test/Er.dat';
[~, name, ~] = fileparts(prof);

%initialize cell array
quant = cell(1, level+1);
%load profiles
raw = load(prof);
x = raw(:, 1);
quant{1} = raw(:, 2);

%numerical gradients
for k = 2:level+1
    quant{k} = gradient(quant{k-1}, x);
end

%use smoothing routine, 0.5% of datapoints, loess method
quant_smooth = smooth2level(quant{1}, x, level, 1e-2, 'loess');

%plot figure
figure('units', 'normalized', 'outerposition', [0, 0, 1, 1]);
axis tight
for k = 1:(level+1)
    subplot(level+1,1,k)
    plot(x, quant{k}, '.r', 'DisplayName', 'original');
    hold on
    plot(x, quant_smooth{k}, '.b', 'DisplayName', 'smoothen');
    hold off
    xlabel('r / cm')
    ylabel('quantity')
    legend('Location', 'best')
    if(k == 1)
        title(name)
    else
        title(['derivative #' num2str(k-1)])
    end
end
set(findall(gcf,'-property','FontSize'), 'FontSize', 18)