function visualizestack(ddir)
% VISUALIZESTACK(ddir)
%
% Plots the events and repeaters for a seismic ocean thermometry
% experiment.
%
% INPUT:
% ddir      directory to the SOT run
% name      name of the SOT run 
%           [default: cindeks(split(ddir, filesep), 'end')]
%
% SEE ALSO:
% QUERYSOT, READQUERYSOT, XCORRSOT
%
% Last modified by spipatprathanporn@ucsd.edu, 10/12/2026

fname = fullfile(ddir, 'staltastack_output.mat');
load(fname, 'station', 'location', 'evs_str')

[lons,lats,elev] = bathymetryglobal;

evs_str.ratio_maxT(isnan(evs_str.ratio_maxT)) = 0;

M = histcounts2(mod(evs_str.PreferredLongitude, 360), ...
    evs_str.PreferredLatitude, 0:1:360, -90:1:90);
M_smooth = conv2(M, ones(3));
M_smooth = M_smooth(2:end-1,2:end-1);

ratio_maxT_raw = zeros(size(M));
for ii = 1:length(evs_str.ratio_maxT)
    i1 = ceil(mod(evs_str.PreferredLongitude(ii),360));
    i2 = ceil(evs_str.PreferredLatitude(ii)+90);
    ratio_maxT_raw(i1,i2) = ratio_maxT_raw(i1,i2) + evs_str.ratio_maxT(ii);
end
ratio_maxT_smooth = conv2(ratio_maxT_raw, ones(3));
ratio_maxT_smooth = ratio_maxT_smooth(2:end-1,2:end-1);

ratio_maxT_avg = ratio_maxT_raw ./ M;
ratio_maxT_avg_smooth = ratio_maxT_smooth ./ M_smooth;

% 1x1 bin-average STA/LTA
figure(1)
clf
set(gcf, 'Units', 'inches', 'Position', [0 1 11 5.5])
ax = subplot('Position', [0.08 0.10 0.88 0.85]);
imagesc(minmax(lons), minmax(lats), elev')
axis xy
axis tight
axis equal
cb = cax2dem([-1 1] * 8000);
delete(cb)
set(ax, 'FontSize', 11, 'TickDir', 'out', 'Box', 'on')
xlabel('longitude')
ylabel('latitude')
xlim([0 360])
ylim([-90 90])
xticks(0:30:360)
yticks(-90:30:90)
title(sprintf('M5.0+ | %s | %d-%d | average STA/LTA | %d events', ...
    station, min(evs_str.dt_origin).Year, max(evs_str.dt_origin).Year, ...
    length(evs_str.ratio_maxT)))

ax2 = addlayeraxes(ax);
im = imagesc(ax2, [0.5 359.5], [-89.5 89.5], ratio_maxT_avg');
set(im, 'AlphaData', ~isnan(ratio_maxT_avg'))
hold on
grid on
plotcont;
h = plotplates;
set(h, 'Color', [0 0 0.75], 'LineWidth', 0.75);
st = scatter(ax2, mod(location(1), 360), location(2), 60, [1 1 0], ...
    'filled', '^', 'MarkerEdgeColor', 'k');
clim([0.5 5])
cbase = [1 1 1; 1 0.8 0.3; 1 0 0; 1 0 1; 0.5 0 1; 0.2 0.8 1];
colormap(ax2, interp1(linspace(0,255,6), cbase, 0:255))
cb = colorbar;
set(cb, 'TickDir', 'out', 'FontSize', 11)
set(get(cb, 'Label'), 'String', 'bin-averaged STA/LTA ratio')
xlim([0 360])
ylim([-90 90])
pause(0.1)
set(ax, 'Position', get(ax2, 'Position'))

set(gcf, 'Renderer', 'painters')
exportgraphics(gcf, fullfile(getenv('EPS'), ...
    sprintf('%s_%s_average_T-phase.pdf', mfilename, ...
    station)), 'ContentType', 'image', 'Resolution', 600, 'Padding', 0);

% smoothed bin-averaged STA/LTA
pause(0.1)
delete(im)
im = imagesc(ax2, [0.5 359.5], [-89.5 89.5], ratio_maxT_avg_smooth');
set(im, 'AlphaData', ~isnan(ratio_maxT_avg_smooth'))
uistack(im, 'bottom')
set(get(cb, 'Label'), 'String', 'smoothed bin-averaged STA/LTA ratio')
exportgraphics(gcf, fullfile(getenv('EPS'), ...
    sprintf('%s_%s_average_T-phase_smoothed.pdf', mfilename, ...
    station)), 'ContentType', 'image', 'Resolution', 600, 'Padding', 0);

% number of Earthquakes
pause(0.1)
delete(im)
im = imagesc(ax2, [0.5 359.5], [-89.5 89.5], M');
set(im, 'AlphaData', (M' > 0))
uistack(im, 'bottom')
clim(ax2, 'auto')
set(get(cb, 'Label'), 'String', 'number of earthquakes in 1x1 degree box')
exportgraphics(gcf, fullfile(getenv('EPS'), ...
    sprintf('%s_%s_earthquake_density.pdf', mfilename, ...
    station)), 'ContentType', 'image', 'Resolution', 600, 'Padding', 0);

% smoothed number of Earthquakes
pause(0.1)
delete(im)
im = imagesc(ax2, [0.5 359.5], [-89.5 89.5], M_smooth');
set(im, 'AlphaData', (M_smooth' > 0))
uistack(im, 'bottom')
clim(ax2, 'auto')
set(get(cb, 'Label'), 'String', 'number of earthquakes in 3x3 degree box')
exportgraphics(gcf, fullfile(getenv('EPS'), ...
    sprintf('%s_%s_earthquake_density_smoothed.pdf', mfilename, ...
    station)), 'ContentType', 'image', 'Resolution', 600, 'Padding', 0);
end