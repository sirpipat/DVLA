function visualizesot(ddir, name)
% VISUALIZESOT(ddir, name)
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
% Last modified by spipatprathanporn@ucsd.edu, 09/19/2026

defval('name', removepath(ddir))

% load relevant variables
fnames = ls2cell(fullfile(ddir, 'querysot_output_*.mat'), 1)';
load(fnames{end}, 'lonlim', 'latlim', 'minmag', 'pstation', 'tstation')
evs_str = readquerysot(ddir);

pkstnm = cindeks(split(pstation, '.'), 2);
tkstnm = cindeks(split(tstation, '.'), 2);
startyear = str2double(indeks(evs_str.PreferredTime{1}, 1:4));
endyear = str2double(indeks(evs_str.PreferredTime{end}, 1:4));

% get the bathymetry map
[lons, lats, elev] = bathymetry([], lonlim, latlim, false);

%% making plots
% all queried events
figure(10)
set(gcf, 'Units', 'inches', 'Position', [0 1 8 8])
ax = subplot('Position', [0.08 0.08 0.84 0.88]);
imagesc(lons, lats, elev')
axis xy
axis tight
axis equal
grid on
cb = cax2dem([-1 1]*8000);
delete(cb);
xlabel('longitude (degrees)')
ylabel('latitude (degrees)')
title(sprintf('M%.1f+ | %d-%d | %d catalog events', minmag, ...
    startyear, endyear, length(evs_str.PreferredTime)))
set(ax, 'Box', 'on', 'TickDir', 'out', 'FontSize', 11)

ax2 = addlayeraxes(ax);
scatter(ax2, mod(evs_str.PreferredLongitude, 360), ...
    evs_str.PreferredLatitude, ...
    3*(evs_str.PreferredMagnitudeValue-2), ...
    years(datetime(evs_str.PreferredTime)-datetime(0,0,0)), 'filled', ...
    'o', 'MarkerEdgeColor', 'k');
colormap(ax2, flipud(kelicol));
cb = colorbar(ax2);
set(cb, 'TickDirection', 'out')
set(get(cb, 'Label'), 'String', 'year')
set(cb, 'FontSize', 11)
set(ax, 'Box', 'on', 'TickDir', 'out', 'FontSize', 11)

set(ax, 'Position', get(ax2, 'Position'))
set(gcf, 'Renderer', 'painters')
figdisp(sprintf('%s_geomap_%s-%s-%s_year_%d-%d', mfilename, name, ...
    pkstnm, tkstnm, startyear, endyear), [], [], 2, [], 'epstopdf')

% repeater graphs
try
    xcorrfname = fullfile(ddir, 'xcorrsot_output.mat');
    load(xcorrfname, 'evs_str', 'dd', 'cc_P', 'cc_T');
    A = and(dd<=60, and(cc_P>=0.9, cc_T>=0.6));
    g = graph(A);

    figure(11)
    set(gcf, 'Units', 'inches', 'Position', [0 1 8 8])
    ax = subplot('Position', [0.08 0.08 0.84 0.88]);
    imagesc(lons, lats, elev')
    axis xy
    axis tight
    axis equal
    grid on
    cb = cax2dem([-1 1]*8000);
    delete(cb);
    xlabel('longitude (degrees)')
    ylabel('latitude (degrees)')
    title(sprintf('Repeater grapth | %d nodes | %d repeaters', ...
        length(evs_str.PreferredTime), height(g.Edges)))
    set(ax, 'Box', 'on', 'TickDir', 'out', 'FontSize', 11)

    ax2 = addlayeraxes(ax);
    for ii = 1:height(g.Edges)
        endnodes = g.Edges(ii,:).EndNodes;
        plot(ax2, evs_str.PreferredLongitude(endnodes), ...
            evs_str.PreferredLatitude(endnodes), 'Color', 'k', ...
            'LineWidth', 0.5)
    end
    hold on
    scatter(evs_str.PreferredLongitude, ...
        evs_str.PreferredLatitude , 6, evs_str.dt_origin.Year, 'filled', ...
        'MarkerEdgeColor', 'k')
    colormap(ax2, jet(endyear-startyear+1))
    clim(ax2, [startyear endyear+1])
    cb = colorbar(ax2);
    set(get(cb, 'Label'), 'String', 'year')
    set(cb, 'FontSize', 11, 'TickDirection', 'both')
    set(ax2, 'Box', 'on', 'TickDir', 'out', 'FontSize', 11)
    if min(get(cb, 'Ticks')) > startyear
        set(cb, 'Ticks', [startyear get(cb, 'Ticks')]);
    end
    if max(get(cb, 'Ticks')) < endyear + 1
        set(cb, 'Ticks', [get(cb, 'Ticks') endyear+1]);
    end

    set(ax, 'Position', get(ax2, 'Position'))
    set(gcf, 'Renderer', 'painters')
    figdisp(sprintf('%s_geomap-graph_%s-%s-%s_year_%d-%d', mfilename, name, ...
        pkstnm, tkstnm, startyear, endyear), [], [], 2, [], 'epstopdf')

    linkaxes([ax ax2]);
catch ME
    ME.getReport
end
end