function targosot_repeaters2022srl
% TARGOSOT_REPEATERS2022SRL
%
% Plots ARGO temperature anomaly along the great-circle paths from
% repeating earthquakes from Yang et al. 2022, SRL, to CTBTO hydroacoustic
% stations.
%
% SEE ALSO:
% TARGOPROFILE, TARGOSOT_SUMATRA
% 
% Last modified by spipatprathanporn@ucsd.edu, 07/31/2026

% repeater table from Yang et al 2022 SRL
D = readtable('~/Documents/IGPP/Research/SRL2022/supplement-mini.csv');
D.Longitude = mod(D.Longitude, 360);
dt = D.Date_ + D.Time;
dt.Format = 'uuuu-MM-dd''T''HH:mm:ss';
% number of repeaters
N = size(D,1) / 2;

% List of hydroacoustic stations run by CTBTO
%         lat    lon is_T-phase
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
hnames = {'CapeLeeuwin'; ...
          'HaidaGwaii'; ...
          'JuanFernandez'; ...
          'Crozet'; ...
          'Guadeloupe'; ...
          'Socorro'; ...
          'Flores'; ...
          'Chagos'; ...
          'TristandaCunha'; ...
          'Ascension'; ...
          'Wake'};

% for each repeaters
for ii = 1:N
    ii_ev1 = 2*ii-1;
    ii_ev2 = 2*ii;

    % skip if any of the event pair is outside 2004-2018
    if D.Date_(ii_ev1).Year < 2004 || D.Date_(ii_ev1).Year > 2018 || ...
            D.Date_(ii_ev2).Year < 2004 || D.Date_(ii_ev2).Year > 2018
        continue
    end

    % Indian basin
    if D.Longitude(ii_ev1) >= 20 && D.Longitude(ii_ev1) < 120
        idx = [1 4 8];
    % West Pacific basin ... Wake Island
    elseif D.Longitude(ii_ev1) >= 120 && mod(D.Longitude(ii_ev1), 360) < 210
        idx = 11;
    % East Pacific basin
    elseif D.Longitude(ii_ev1) >= 210 && D.Longitude(ii_ev1) < 300
        % Mexico
        if D.Latitude(ii_ev1) >= -10
            idx = 6;
        % Chile
        else
            idx = 3;
        end
    % South Atlantic
    else
        idx = [9 10];
    end

    stlo = hloc(idx, 2);
    stla = hloc(idx, 1);
    kstnm = hnames(idx);

    for jj = 1:length(stlo)
        % Event 1
        [x,z,T,dT] = targoprofile([D.Longitude(ii_ev1) D.Latitude(ii_ev1)], ...
                [stlo(jj) stla(jj)], 120, dt(ii_ev1));
        figure(10)
        clf
        set(gcf, 'Units', 'inches', 'Position', [0 1 9 8])
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
            D.Latitude(ii_ev1), D.Longitude(ii_ev1), stla(jj), stlo(jj), dt(ii_ev1)))
        set(get(cb, 'Label'), 'String', 'ARGO mean temperature ({\circ}C)')
        set(gca, 'FontSize', 11, 'Box', 'on', 'TickDir', 'out')

        subplot(212)
        dT_log = log10(abs(dT));
        dT_sgn = sign(dT);
        dT_plot = max(dT_log + 2, 0) .* dT_sgn;
        im = imagesc(x, z, dT_plot);
        set(im, 'AlphaData', ~isnan(dT_plot))
        colormap(flipud(kelicol))
        cl = clim;
        %clim([-1 1] * max(cl))
        clim([-3 3])
        cb = colorbar;
        grid on
        xlabel('distance along the great circle path (km)')
        ylabel('depth (m)')
        title('ARGO temperature anomaly')
        subtitle(sprintf(['Left: (lat,lon) = (%.2f,%.2f) | Right: ' ...
            '(lat,lon) = (%.2f,%.2f) | Event: %s'], ...
            D.Latitude(ii_ev1), D.Longitude(ii_ev1), stla(jj), stlo(jj), dt(ii_ev1)))
        set(get(cb, 'Label'), 'String', 'ARGO temperature anomaly ({\circ}C)')
        set(cb, 'TickLabels', {-10, -1, -0.1, '|dT| < 0.01', 0.1, 1, 10})
        set(gca, 'FontSize', 11, 'Box', 'on', 'TickDir', 'out')

        set(gcf, 'Renderer', 'painters')
        figdisp(sprintf('%s_N%03d-1_%s_%s', mfilename, ii, ...
            replace(string(dt(ii_ev1)), ':', '-'), ...
            kstnm{jj}), [], [], 2, [], 'epstopdf')

        ev1.x = x; ev1.z = z; ev1.dT = dT;

        % Event 2
        [x,z,T,dT] = targoprofile([D.Longitude(ii_ev2) D.Latitude(ii_ev2)], ...
                [stlo(jj) stla(jj)], 120, dt(ii_ev2));
        figure(10)
        clf
        set(gcf, 'Units', 'inches', 'Position', [0 1 9 8])
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
            D.Latitude(ii_ev2), D.Longitude(ii_ev2), stla(jj), stlo(jj), dt(ii_ev2)))
        set(get(cb, 'Label'), 'String', 'ARGO mean temperature ({\circ}C)')
        set(gca, 'FontSize', 11, 'Box', 'on', 'TickDir', 'out')

        subplot(212)
        dT_log = log10(abs(dT));
        dT_sgn = sign(dT);
        dT_plot = max(dT_log + 2, 0) .* dT_sgn;
        im = imagesc(x, z, dT_plot);
        set(im, 'AlphaData', ~isnan(dT_plot))
        colormap(flipud(kelicol))
        cl = clim;
        %clim([-1 1] * max(cl))
        clim([-3 3])
        cb = colorbar;
        grid on
        xlabel('distance along the great circle path (km)')
        ylabel('depth (m)')
        title('ARGO temperature anomaly')
        subtitle(sprintf(['Left: (lat,lon) = (%.2f,%.2f) | Right: ' ...
            '(lat,lon) = (%.2f,%.2f) | Event: %s'], ...
            D.Latitude(ii_ev2), D.Longitude(ii_ev2), stla(jj), stlo(jj), dt(ii_ev2)))
        set(get(cb, 'Label'), 'String', 'ARGO temperature anomaly ({\circ}C)')
        set(cb, 'TickLabels', {-10, -1, -0.1, '|dT| < 0.01', 0.1, 1, 10})
        set(gca, 'FontSize', 11, 'Box', 'on', 'TickDir', 'out')

        set(gcf, 'Renderer', 'painters')
        figdisp(sprintf('%s_N%03d-2_%s_%s', mfilename, ii, ...
            replace(string(dt(ii_ev2)), ':', '-'), ...
            kstnm{jj}), [], [], 2, [], 'epstopdf')

        ev2.x = x; ev2.z = z; ev2.dT = dT;
        evs = [ev1 ev2];

        % Path coherence
        for kk = 1:2
            evs(kk).dT_mean = mean(evs(kk).dT, 2, 'omitmissing');
            evs(kk).dT_std  = std(evs(kk).dT,[],2,'omitmissing');
            evs(kk).dT_z = (evs(kk).dT - evs(kk).dT_mean) ./ evs(kk).dT_std;
        end
        dT_cc = evs(1).dT_z .* evs(2).dT_z;
        dT_cc = mean(dT_cc,2,'omitmissing');

        figure(12)
        clf
        set(gcf, 'Units', 'inches', 'Position', [0 1 9 8])
        subplot('Position', [0.08 0.94 0.84 0.01])
        title(sprintf('N:%03d | %s -- %s | %s | %d km', ii, dt(ii_ev1), ...
            dt(ii_ev2), kstnm{jj}, round(evs(2).x(end))))
        nolabels(gca, 3)
        set(gca, 'Color', 'none', 'FontSize', 12)
        set(get(gca, 'XAxis'), 'Visible', 'off')
        set(get(gca, 'YAxis'), 'Visible', 'off')

        subplot(131)
        plot(evs(1).dT_mean, evs(1).z, 'LineWidth', 1, 'Color', 'b')
        hold on
        plot(evs(2).dT_mean, evs(2).z, 'LineWidth', 1, 'Color', 'r')
        axis ij
        grid on
        legend('event 1', 'event 2', 'Location', 'southeast')
        xlabel('average along-path anomaly (C)')
        ylabel('depth (m)')
        title('Average along-path anomaly')
        set(gca, 'Box', 'on', 'TickDir', 'out', 'FontSize', 10)

        subplot(132)
        plot(evs(2).dT_mean - evs(1).dT_mean, evs(2).z, 'LineWidth', 1, 'Color', 'k')
        axis ij
        grid on
        xlim([-1 1] * 0.2)
        ylim([0 2000])
        xticks((-1:0.5:1) * 0.2)
        xlabel('average along-path anomaly difference (C)')
        ylabel('depth (m)')
        title('Average along-path anomaly difference')
        set(gca, 'Box', 'on', 'TickDir', 'out', 'FontSize', 10)

        subplot(133)
        plot(dT_cc, evs(2).z, 'LineWidth', 1, 'Color', 'k')
        axis ij
        grid on
        xlim([-1 1])
        xticks(-1:0.5:1)
        xlabel('correlation coefficient')
        ylabel('depth (m)')
        title('Coherence')
        set(gca, 'Box', 'on', 'TickDir', 'out', 'FontSize', 10)

        set(gcf, 'Renderer', 'painters')
        figdisp(sprintf('%s_N%03d_coherence_%s', mfilename, ii, ...
            kstnm{jj}), [], [], 2, [], 'epstopdf')
    end
end
% MORE IDEAS:
% 1. Event time difference vs coherence or average along-path dT difference
% 2. Same as 1 but only event tiem difference within a year
end