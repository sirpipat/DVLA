function repeaterpaths
% REPEATERPATHS
%
% Plots CTBTO hydroacoustic stations and repeaters locations from various
% compilations
% 
% Last modified by spipatprathanporn@ucsd.edu, 09/08/2026

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

% which station is a hydrophone station
wh = (hloc(:,3) == 0);
% colors of the station type
cmap = [1 0.75 0.15; 0 0.75 0.2];


% repeater table from Yang et al 2022 SRL
T = readtable('~/Documents/IGPP/Research/SRL2022/supplement-mini.csv');
wht = and(T.Date.Year >= 1900, T.Date.Year <= 2090);
T.Longitude = mod(T.Longitude, 360);

% repeater table from 

% plots the world map
figure(1)
set(gcf, 'Units', 'inches', 'Position', [0 1 8 5])
clf
[~, handl_cont] = plotcont;
set(handl_cont, 'HandleVisibility', 'off')

hold on
handl_plate = plotplates;
set(handl_plate, 'Color', 'r', 'HandleVisibility', 'off')

% plots the repeater epicenters
rp_2022srl = scatter(T.Longitude(wht), T.Latitude(wht), 40, ...
    T.Date.Year(wht), 'filled', 'p', 'MarkerEdgeColor', 'k', 'LineWidth', 0.25);
% plots hydrophone station locations
hy = scatter(mod(hloc(wh,2), 360), hloc(wh,1), 60, cmap(1,:), ...
    'filled', '^', 'MarkerEdgeColor', 'k');
% plots T-wave station locations
tw = scatter(mod(hloc(~wh,2), 360), hloc(~wh,1), 60, cmap(2,:), ...
    'filled', 'v', 'MarkerEdgeColor', 'k');

cb = colorbar;
set(cb, 'TickDir', 'out')
set(get(cb, 'Label'), 'String', 'Year')

% plot paths
for ii = 1:length(T.Longitude)/2
    ii_ev1 = 2*ii-1;
    ii_ev2 = 2*ii;
    % skip if any of the event pair is outside 2004-2018
    if ~wht(ii_ev1) || ~wht(ii_ev2)
        continue
    end

    idx = getimsindex(T.Longitude(ii_ev1), T.Latitude(ii_ev1));
    
    stlo = hloc(idx, 2);
    stla = hloc(idx, 1);

    for jj = 1:length(stlo)
        plottrack(gca, [T.Longitude(ii_ev1) T.Latitude(ii_ev1)], [stlo(jj) stla(jj)], 0, 120, 'Color', [0.5 0.5 1], 'LineWidth', 0.25);
    end
end

legend('Yang et al. 2022 SRL', 'hydrophone station', ...
    'T-phase station (island)', 'Location', 'southoutside')
axis xy
axis tight
axis equal
grid on

uistack(handl_cont, 'top')
uistack(handl_plate, 'top')
uistack(rp_2022srl, 'top')
uistack(hy, 'top')
uistack(tw, 'top')

xlabel('longitude (degree)')
ylabel('latitude (degree)')
xticks(0:30:360)
colormap("jet")
title('Repeater paths')
set(gca, 'TickDir', 'out', 'Box', 'on', 'FontSize', 12)
set(gcf, 'Renderer', 'painters')
figdisp(mfilename, [], [], 2, [], 'epstopdf')
end

function idx = getimsindex(lon, lat)
% Indian basin
if lon >= 20 && lon < 120
    idx = [1 4 8];
% West Pacific basin ... Wake Island
elseif lon >= 120 && mod(lon, 360) < 210
    idx = 11;
% East Pacific basin
elseif lon >= 210 && lon < 300
    % Mexico
    if lat >= -10
        idx = 6;
    % Chile
    else
        idx = 3;
    end
% South Atlantic
else
    idx = [9 10];
end
end