function ax2 = addlayeraxes(ax)
% ax2 = ADDLAYERAXES(ax)
%
% Adds a transparent axes on top of the original axes, making a new layer 
% on top. It is good for plotting two or more objects with different
% colormap in the same axes.
%
% INPUT:
% ax        base axes
%
% OUTPUT:
% ax2       axes of the new layer
%
% SEE ALSO:
% DOUBLEAXES
%
% Last modified by spipatprathanporn@ucsd.edu, 09/08/2026

ax2 = axes('Parent', get(ax, 'Parent'), ...
    'Position', get(ax, 'Position'), ...
    'XLimMode', 'manual', ...
    'XLim', get(ax, 'XLim'), ...
    'YLimMode', 'manual', ...
    'YLim', get(ax, 'YLim'), ...
    'DataAspectRatioMode', 'manual', ...
    'XTick', [], 'YTick', [], 'Color', 'none');
nolabels(ax2, 3);
linkaxes([ax ax2]);
hold on
end
