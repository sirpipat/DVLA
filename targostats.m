function targostats
% TARGOSTATS
%
% Plots standard deviation of ARGO temperature anomaly at any given depth
% with CTBTO hydroacoustic stations and repeaters locations from Yang et
% al. 2022 SRL
% 
% Last modified by spipatprathanporn@ucsd.edu, 07/31/2026

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
T = readtable('~/Documents/IGPP/Research/SRL2021/supplement-mini.csv');
wht = and(T.Date_.Year >= 2004, T.Date_.Year <= 2018);

% reads ARGO gridded temperature
ncfile = fullfile(getenv('IFILES'), 'EARTHMODELS', 'PHYSICAL', ...
    'ARGO', 'RG_ArgoClim_Temperature_2019.nc');

% reads dimensions of the grid
LONGITUDE = ncread(ncfile, 'LONGITUDE');
LATITUDE = ncread(ncfile, 'LATITUDE');
PRESSURE = ncread(ncfile, 'PRESSURE');
TIME = ncread(ncfile, 'TIME');

nlon = length(LONGITUDE);
nlat = length(LATITUDE);
npres = length(PRESSURE);
ntime = length(TIME);

% loops over pressure (depth)
for ii = 1:npres
    % plots the ocean map of STD of temperature anomaly
    figure(1)
    set(gcf, 'Units', 'inches', 'Position', [0 1 8 5])
    clf
    ARGO_TEMPERATURE_ANOMALY = ncread(ncfile, ...
        'ARGO_TEMPERATURE_ANOMALY', [1 1 ii 1], [Inf Inf 1 Inf]);
    ARGO_TEMPERATURE_STD = std(ARGO_TEMPERATURE_ANOMALY, [], 4, ...
        "omitmissing");
    im = imagesc(LONGITUDE, LATITUDE, ARGO_TEMPERATURE_STD');
    set(im, 'AlphaData', ~isnan(ARGO_TEMPERATURE_STD'))
    hold on
    % plots the repeater epicenters
    scatter(mod(T.Longitude(wht)-20,360)+20, T.Latitude(wht), 30, ...
        [0 0 0], 'filled', 'p', 'MarkerEdgeColor', 'w', 'LineWidth', 0.25)
    % plots hydrophone station locations
    scatter(mod(hloc(wh,2), 360), hloc(wh,1), 60, cmap(1,:), ...
        'filled', '^', 'MarkerEdgeColor', 'k')
    % plots T-wave station locations
    scatter(mod(hloc(~wh,2), 360), hloc(~wh,1), 60, cmap(2,:), ...
        'filled', 'v', 'MarkerEdgeColor', 'k')
    legend('repeaters', 'hydrophone station', ...
        'T-phase station (island)', 'Location', 'southoutside')
    axis xy
    axis tight
    axis equal
    grid on
    xlabel('longitude (degree)')
    ylabel('latitude (degree)')
    xticks(0:30:360)
    cb = colorbar;
    set(get(cb, 'Label'), 'String', 'Temperature anomaly std ({^\circ}C)')
    set(cb, 'TickDirection', 'out')
    colormap("jet")
    title(sprintf('ARGO temperature anomaly at %d mbar (%g m)', ...
        PRESSURE(ii)*100, PRESSURE(ii)))
    set(gca, 'TickDir', 'out', 'Box', 'on', 'FontSize', 12)
    set(gcf, 'Renderer', 'painters')
    exportgraphics(gcf, fullfile(getenv('EPS'), ...
        sprintf('%s_std_%06d-mbar.pdf', mfilename, PRESSURE(ii)*100)), ...
        'ContentType', 'image', 'Resolution', 600, 'Padding', 0)
    % figdisp(sprintf('%s_std_%06d-mbar', mfilename, PRESSURE(ii)*100), [], ...
    %     [], 2, [], 'epstopdf')

    % plots the ocean map of temperature anomaly trend
    figure(2)
    set(gcf, 'Units', 'inches', 'Position', [0 1 8 5])
    clf
    % linear trend fitting
    % time in years
    t = ((1:ntime)' - (1+ntime)/2 ) / 12;
    G = [t ones(size(t))];
    % combine LONGITUDE and LATITUDE to one dimension
    d = reshape(squeeze(ARGO_TEMPERATURE_ANOMALY), nlon*nlat, ntime)';
    m = G \ d;
    ARGO_TEMPERATURE_TREND = reshape(m(1,:), nlon, nlat);
    im = imagesc(LONGITUDE, LATITUDE, ARGO_TEMPERATURE_TREND');
    set(im, 'AlphaData', ~isnan(ARGO_TEMPERATURE_TREND'))
    hold on
    % plots the repeater epicenters
    scatter(mod(T.Longitude(wht)-20,360)+20, T.Latitude(wht), 30, ...
        [0 0 0], filled', 'p')
    % plots hydrophone station locations
    scatter(mod(hloc(wh,2), 360), hloc(wh,1), 60, cmap(1,:), ...
        'filled', '^', 'MarkerEdgeColor', 'k')
    % plots T-wave station locations
    scatter(mod(hloc(~wh,2), 360), hloc(~wh,1), 60, cmap(2,:), ...
        'filled', 'v', 'MarkerEdgeColor', 'k')
    legend('repeaters', 'hydrophone station', ...
        'T-phase station (island)', 'Location', 'southoutside')
    axis xy
    axis tight
    axis equal
    grid on
    xlabel('longitude (degree)')
    ylabel('latitude (degree)')
    xticks(0:30:360)
    cb = colorbar;
    set(get(cb, 'Label'), 'String', ['Temperature anomaly trend ' ...
        '({^\circ}C/yr)'])
    set(cb, 'TickDirection', 'out')
    colormap("jet")
    cl = clim;
    clim([-1 1] * max(cl))
    title(sprintf('ARGO temperature anomaly trend at %d mbar (%g m)', ...
        PRESSURE(ii)*100, PRESSURE(ii)))
    set(gca, 'TickDir', 'out', 'Box', 'on', 'FontSize', 12)
    set(gcf, 'Renderer', 'painters')
    exportgraphics(gcf, fullfile(getenv('EPS'), ...
        sprintf('%s_trend_%06d-mbar.pdf', mfilename, ...
        PRESSURE(ii)*100)), 'ContentType', 'image', 'Resolution', 600, ...
        'Padding', 0)
    % figdisp(sprintf('%s_trend_%06d-mbar', mfilename, PRESSURE(ii)*100), [], ...
    %     [], 2, [], 'epstopdf')
end
end