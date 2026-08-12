function imshydrophoneloc
% IMSHYDROPHONELOC
%
% Plots the location of IMS hydroacoustic stations on a world map.
%
% Last modified by spipatprathanporn@ucsd.edu, 08/12/2026

% lat lon is_T-phase
hloc = [-34.3  115.2 0; ...   % Cape Leeuwin, Australia
         53.3 -132.5 1; ...   % Haida Gwaii, Canada
        -33.6  -78.8 0; ...   % Juan Fernandez Island, Chile
        -46.4   51.9 0; ...   % Crozet Islands, France
         16.3  -61.1 1; ...   % Guadeloupe, France
         18.7 -110.9 1; ...   % Socorro Island, Mexico
         39.4  -31.2 1; ...   % Flores, Portugal
         -7.3   72.4 0; ...   % BIOT/Chagos Archipelago, UK
        -37.1  -12.3 1; ...   % Tristan da Cunha, UK
         -8.0  -14.4 0; ...   % Ascension, UK
         19.3  166.6 0];      % Wake Island, US

% colors of the hydroacoustic stations
cmap = [1 0.75 0.15; 0 0.75 0];

% makes a figure
figure(9)
clf
set(gcf, 'Units', 'inches', 'Position', [0 1 8 6])
hold on
% plots the constlines
[~, co] = plotcont;
set(co, 'HandleVisibility', 'off')
% plots the plate boundaries
pl = plotplates;
set(pl, 'Color', 'r')
set(pl, 'HandleVisibility', 'off')
% plots the IMS hydroacustic stations
wh = (hloc(:,3)==0);
scatter(mod(hloc(wh,2), 360), hloc(wh,1), 60, cmap(1,:), 'filled', ...
    '^', 'MarkerEdgeColor', 'k')
scatter(mod(hloc(~wh,2), 360), hloc(~wh,1), 60, cmap(2,:), 'filled', ...
    'v', 'MarkerEdgeColor', 'k')
legend('hydrophone station', 'T-phase station (island)')
% makes the plot look pretty
grid on
xlim([0 360])
ylim([-90 90])
xticks(0:30:360)
yticks(-90:30:90)
xlabel('longitude')
ylabel('latitude')
title('IMS hydroacoustic stations')
set(gca, 'Box', 'on', 'TickDir', 'out', 'FontSize', 12)
set(gcf, 'Renderer', 'painters')
% saves the figure
figdisp(mfilename, [], [], 2, [], 'epstopdf')
end