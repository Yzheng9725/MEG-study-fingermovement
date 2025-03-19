function wxyz_decorFigure(hFig, varargin)
% WXYZ_DECORFIGURE This function implements a figure that creates a tiledlayout
% containing the specified number of grids. 
% 
% This function takes a handle to figure and the corresponding arguments to
% modify the figure.
% 
% hFig  - the first input parameter, should be a positive integer.
% 
% Other arguments should come in key-value pairs and can include
% 'FontSize', 'FontName', 'FontWeight', 'FontScale', 'BoxWidth', 'Box'.
%
% example:
%   wxyz_decorFigure(hFig, 'FontSize', 12, 'Box', 'off');
% 
% Author: wxyz
% Version: 1.0
% Last revision date : 2024-01-05

% Code implementation section
% Default value
params.FontSize = 12;
params.FontName = 'Calibri';
params.FontWeight = 'bold';
params.FontScale = 1.4;
params.BoxWidth = 1;
params.Box = 'on';

for i = 1:2:length(varargin) % parse varargin
    key = varargin{i};
    if isfield(params, key)
        params.(key) = varargin{i+1};
    else
        error(['Unrecognized key: ', key]);
    end
%     value = varargin{i+1};
%     switch key
%         case 'FontSize'
%             FontSize = value;
%         case 'FontName'
%             FontName = value;
%         case 'FontWeight'
%             FontWeight = value;
%         case 'FontScale'
%             FontScale = value;
%         case 'Box'
%             Box = value;
%         case 'BoxWidth'
%             BoxWidth = value;
%         otherwise
%             error(['Unrecognized key: ', key]);
%     end
end

%% Font
set(hFig, 'Color', 'w');
set(findobj(hFig, 'Type', 'Axes'), ...
    'FontName', params.FontName,...
    'FontSize', params.FontSize,...
    'FontWeight', params.FontWeight,...
    'TitleFontWeight', params.FontWeight, ...
    'LabelFontSizeMultiplier', params.FontScale,...
    'TitleFontSizeMultiplier', params.FontScale);

if ~isempty(findobj(hFig, 'Type', 'Legend'))
    set(findobj(hFig, 'Type', 'Legend'), 'FontSize', params.FontSize*1, 'FontWeight', params.FontWeight);
end

if ~isempty(findobj(hFig, 'Type', 'Text'))
    set(findobj(hFig, 'Type', 'Text'), 'FontName', params.FontName, 'FontSize', params.FontSize*1);
end

if ~isempty(findobj(hFig, 'Type', 'ColorBar'))
    set(findobj(hFig, 'Type', 'ColorBar'), 'FontSize', params.FontSize*0.9);
    cb = findall(hFig, 'Type', 'Colorbar');
    for c = 1:numel(cb)
        set(cb(c).Label, 'FontName', params.FontName, 'FontSize', params.FontSize, 'FontWeight', params.FontWeight, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'cap');
    end
end

if ~isempty(findobj(hFig, 'Type', 'tiledlayout'))
    set(findobj(hFig, 'Type', 'tiledlayout').Title, 'FontSize', params.FontSize*params.FontScale,...
        'FontName', params.FontName, 'FontWeight', params.FontWeight);
end
if ~isempty(findall(hFig, 'Type', 'textbox')) % Annotation
    set(findall(hFig, 'Type', 'textbox'), 'FontSize', params.FontSize*params.FontScale,...
        'FontName', params.FontName, 'FontWeight', params.FontWeight, 'FitBoxToText', true, 'EdgeColor', 'none');
end

%% Box
set(findobj(hFig, 'Type', 'Axes'), 'Box', params.Box); % axes box
set(findobj(hFig, 'Type', 'Axes'), 'LineWidth', params.BoxWidth); % axes box line


