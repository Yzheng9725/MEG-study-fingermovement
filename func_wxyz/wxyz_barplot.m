function [hbar, hlgd] = wxyz_barplot(data, err, color, varargin) % dat, err, lgd, xaxistick, yaxislabel, gcatitle, color, isPlot0err)
% WXYZ_BARPLOT This function implements the ability to plot bar graphs and
% error bars.
% 
% This function takes data, err and color, as well as other arguments, and 
% then returns handles to the drawn bars and legend.
% 
% data  - used to plot the bars should be an array of m*n, where m is the
%         group and n is the number of conditions in each group. 
% err   - used to plot the error bars has the same dimensions as data.
% color - used for plotting, which should be a matrix of n*3.
%
% Other arguments should come in key-value pairs and can include 'legend',
% 'xtick', 'xlabel', 'ylabel', 'title', 'show0err', 'BarWidth',
% 'errBarStyle', 'errBarWidth', 'errBarColor'.
% 
% hbar  - the created bar plot handle.
% hlgd  - the created elgend handle.
%
% example:
%   [hbar, hlgd] = wxyz_barplot(data, err, color, 'xtick', {'G1', 'G2', 'G3'});
% 
% Author: wxyz
% Version: 1.0
% Last revision date : 2024-01-08

% dat: m*n. m->group, n->barnum per group.

% Code implementation section
% Everything is added to the current figure
holdflag = ishold;
if ~holdflag
  hold on
end

% Default value
params.data         = data;
params.err          = err;
params.color        = color;
params.legend       = '';
params.showlegend   = true;
params.xtick        = '';
params.xlabel       = '';
params.ylabel       = '';
params.title        = '';
params.show0err     = true;
params.BarWidth     = 1;
params.errBarStyle  = '-';
params.errBarWidth  = 1;
params.errBarColor  = 'k';
params.facealpha    = 1;

if ~isempty(varargin)
    for i = 1:2:length(varargin) % parse varargin
        key = varargin{i};
        if isfield(params, key)
            params.(key) = varargin{i+1};
        else
            error(['Unrecognized key: ', key]);
        end
    end
end

if size(params.color, 1) ~= size(params.data, 2)
    error('Color number should be euqal with bar number.');
end

hbar = bar(params.data, 'BarWidth', params.BarWidth);
for i = 1:size(params.data, 2)
    hbar(i).FaceColor = params.color(i, :);
    hbar(i).FaceAlpha = params.facealpha;
    if ~isempty(params.err)
        if params.show0err
            errorbar(hbar(i).XEndPoints, params.data(:, i), params.err(:, i), 'LineStyle', 'none', 'Color', params.errBarColor, 'LineWidth', params.errBarWidth);
        else
            if ~isempty(find(params.err(:,i)~=0, 1))
                errorbar(hbar(i).XEndPoints(params.err(:,i)~=0), data(params.err(:,i)~=0, i), params.err(params.err(:,i)~=0, i), 'LineStyle', params.errBarStyle, 'Color', params.errBarColor, 'LineWidth', params.errBarWidth);
            end
        end
    end
end
if params.showlegend
    if isempty(params.legend)
        params.legend = [];
        for i = 1:size(params.data, 2)
            params.legend{i} = strcat('Condition', num2str(i));
        end
    end
    hlgd = legend(hbar, params.legend);
else
    hlgd = [];
end

if isempty(params.xtick)
    params.xtick = [];
    for i = 1:size(params.data, 1)
        params.xtick{i} = strcat('Group', num2str(i));
    end
end
set(gca, 'XTick', 1:size(params.data, 1), 'XTickLabel', params.xtick);

xlabel(params.xlabel);
ylabel(params.ylabel);
title(params.title);

if ~holdflag
  hold off
end

if ~nargout
  clear hbar hlgd
end