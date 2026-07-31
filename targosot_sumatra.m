function targosot_sumatra
% TARGOSOT_SUMATRA
%
% Plots ARGO temperature anomaly along the great-circle paths from
% repeating earthquakes in Sumatra to a seismic station in Diego Garcia
% island.
%
% SEE ALSO:
% TARGOPROFILE, TARGOSOT_REPEATERS2022SRL
%
% Last modified by spipatprathanporn@ucsd.edu, 07/31/2026

% loads computed P- and T-wave cross-correlation as well as event
% parameters from XCORRSOT
fname = fullfile('/Users/spipatprathanporn/research/IFILES/MATFILES', ...
    'sot_sumatra-psi-dgar_20260509.mat');
load(fname, 'cc_P', 'cc_T', 'dd', 'evs_str' );
A = and(and(cc_T >= 0.6, cc_P >= 0.9), dd <= 60);

% identifies events of the repeaters
wh = any(A,1);
dt_used = datetime(evs_str.PreferredTime(wh));
evlo_used = evs_str.PreferredLongitude(wh);
evla_used = evs_str.PreferredLatitude(wh);
id_used = evs_str.PublicId(wh);

% Diego Garcia location
stlo = 72.4525;
stla = -7.4121;

% loops over events
for ii = 1:length(evlo_used)
    % Makes a cross-section of ARGO mean temperature and temperature
    % anomaly
    [x,z,T,dT] = targoprofile([evlo_used(ii) evla_used(ii)], ...
        [stlo stla], 120, dt_used(ii));

    % plot
    figure(10)
    clf
    set(gcf, 'Units', 'inches', 'Position', [0 1 9 8])

    % plot the mean temperature
    subplot(211)
    im = imagesc(x, z, T);
    set(im, 'AlphaData', ~isnan(T))
    colormap(flipud(kelicol))
    clim([0 30])
    cb = colorbar;
    grid on
    xlabel('distance along the great circle path (km)')
    ylabel('depth (m)')
    title('ARGO mean temperature')
    subtitle(sprintf(['Left: (lat,lon) = (%.2f,%.2f) | Right: ' ...
        '(lat,lon) = (%.2f,%.2f) | Event: %s'], ...
        evla_used(ii), evlo_used(ii), stla, stlo, dt_used(ii)))
    set(get(cb, 'Label'), 'String', 'ARGO mean temperature ({\circ}C)')
    set(gca, 'FontSize', 11, 'Box', 'on', 'TickDir', 'out')

    % plot the temperature anomaly in a signed log scale with a cut off at
    % 0.01 degree celcius
    subplot(212)
    dT_log = log10(abs(dT));
    dT_sgn = sign(dT);
    dT_plot = max(dT_log + 2, 0) .* dT_sgn;
    im = imagesc(x, z, dT_plot);
    set(im, 'AlphaData', ~isnan(dT_plot))
    colormap(flipud(kelicol))
    clim([-3 3])
    cb = colorbar;
    grid on
    xlabel('distance along the great circle path (km)')
    ylabel('depth (m)')
    title('ARGO temperature anomaly')
    subtitle(sprintf(['Left: (lat,lon) = (%.2f,%.2f) | Right: ' ...
        '(lat,lon) = (%.2f,%.2f) | Event: %s'], ...
        evla_used(ii), evlo_used(ii), stla, stlo, dt_used(ii)))
    set(get(cb, 'Label'), 'String', 'ARGO temperature anomaly ({\circ}C)')
    set(cb, 'TickLabels', {-10, -1, -0.1, '|dT| < 0.01', 0.1, 1, 10})
    set(gca, 'FontSize', 11, 'Box', 'on', 'TickDir', 'out')

    set(gcf, 'Renderer', 'painters')

    % save the figure
    figdisp(sprintf('%s_%s', mfilename, ...
        cindeks(split(id_used{ii}, '='), 'end')), [], [], 2, [], ...
        'epstopdf')
end
end