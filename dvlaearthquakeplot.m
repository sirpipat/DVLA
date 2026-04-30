function dvlaearthquakeplot(ev)
% DVLAEARTHQUAKEPLOT(ev)
%
% Plots selected hydrophone sections at the hydrophone moorings with the
% distance and predicted P- and T- phase arrival time. This function is for
% figure generation only.
%
% INPUT:
% ev        deployment day number of the earthquake origin time:
%           options:
%           + 373
%           + 291
%           + 294 - plot the polarmap instead of seismograms
%
% Last modified by spipatprathanporn@ucsd.edu, 04/30/2026

% sampling rate
fs = 20;
% resp removal output
resp_output_fmt = 'velocity';

if ev == 373
    ev373 = irisFetch.Events('MinimumMagnitude', 4.5, ...
        'latitude', 86.228, 'longitude', 35.748, 'maximumRadius', 1, ...
        'startTime', '2020-01-08T12:00:00', ...
        'endTime', '2020-01-08T12:59:59')';
    dt_origin = datetime(ev373.PreferredTime, ...
        'Format', 'dd-MMM-uuuu HH:mm:ss.SSSSSS', 'TimeZone', 'UTC');
    evlo = ev373.PreferredLongitude;
    evla = ev373.PreferredLatitude;
    evdp = ev373.PreferredDepth;
    mag = ev373.PreferredMagnitudeValue;
    
    % seismogram plot
    fig3 = figure(3);
    clf
    set(fig3, 'Units', 'inches', 'Position', [0 1 6 8])
    subplot('Position', [0.12 0.08 0.84 0.82])
    hold on

    % Land stations
    ddir = '/Users/spipatprathanporn/research/IFILES/SEISMOQUERY/Event_11167897/';

    station_list = {'II.ALE', 'TA.A21K', 'TA.C26K', ...
        'TA.D28M', 'TA.A36M'};
    
    for ii = 1:length(station_list)
        station = station_list{ii};
        search_str1 = sprintf('%s*BHE*.sac', station);
        search_str2 = sprintf('%s*BH2*.sac', station);
        try
            fname = cindeks(ls2cell(fullfile(ddir, search_str1), 1), 1);
        catch
            fname = cindeks(ls2cell(fullfile(ddir, search_str2), 1), 1);
        end
        [x, HdrData] = readsac(fname);
        [~, ~, ~, fs, ~, dt, ~] = gethdrinfo(HdrData);
        stlo = HdrData.STLO;
        stla = HdrData.STLA;
    
        sacpzfile = [fname, 'pz'];
        x = detrend(x, 2) .* shanning(length(x), 0.016);
        x = real(transfer(x, 1/fs, [0.01 0.02 8 10], resp_output_fmt, sacpzfile, 'sacpz'));
        x = bandpass(x, fs, 2, 6, 2, 1, 'butter', 'linear');
        t = seconds(dt - dt_origin);
        [distkm, distdeg] = grcdist([evlo evla], [stlo stla]);
        if distdeg < 20
            plot(t, x / max(abs(x)) * 2 + distdeg, 'Color', [0.65 0.65 0.65], 'LineWidth', 1)
        else
            plot(t, x / max(abs(x)) / 2 + distdeg, 'Color', [0.65 0.65 0.65], 'LineWidth', 1)
        end
        text(2000, distdeg+0.35, station, 'FontSize', 9, ...
            'Color', [0.3 0.3 0.3])
    end

    % NERSC1
    info_mooring = mooringinfo('NERSC1');
    stlo = info_mooring.rreflon;
    stla = info_mooring.rreflat;
    [~, wh] = min(abs(info_mooring.hmdep-1000));
    [t, x, jdn, num, depth] = getdvlaseis('NERSC1', info_mooring.snx(wh), ...
        '2020-01-08 12:07:00', '2020-01-08 12:37:00', fs);
    dt = datetime(jdn, 'ConvertFrom', 'datenum', ...
        'Format', 'uuuu-MM-dd''T''HH:mm:ss.SSSSSS', ...
        'TimeZone', 'UTC') + seconds(t);
    t_nersc1 = seconds(dt - dt_origin);
    x = bandpass(x, fs, 0.5, 5, 2, 1, 'butter', 'linear');
    [distkm, distdeg] = grcdist([evlo evla], [stlo stla]);
    plot(t_nersc1, x / max(abs(x)) + distdeg, ...
        'Color', 'k', 'LineWidth', 1)
    scatter(1800, distdeg+0.15, 40, [1 0 0], 'filled', 'v', 'MarkerEdgeColor', 'k');
    text(1850, distdeg+0.1, 'NERSC1', 'Color', 'k', 'FontSize', 9)

    % NERSC2
    info_mooring = mooringinfo('NERSC2');
    stlo = info_mooring.rreflon;
    stla = info_mooring.rreflat;
    [~, wh] = min(abs(info_mooring.hmdep-1000));
    [t, x, jdn, num, depth] = getdvlaseis('NERSC2', info_mooring.snx(wh), ...
        '2020-01-08 12:07:00', '2020-01-08 12:37:00', fs);
    dt = datetime(jdn, 'ConvertFrom', 'datenum', ...
        'Format', 'uuuu-MM-dd''T''HH:mm:ss.SSSSSS', ...
        'TimeZone', 'UTC') + seconds(t);
    t_nersc1 = seconds(dt - dt_origin);
    x = bandpass(x, fs, 0.5, 5, 2, 1, 'butter', 'linear');
    [distkm, distdeg] = grcdist([evlo evla], [stlo stla]);
    plot(t_nersc1, x / max(abs(x)) + distdeg, ...
        'Color', 'k', 'LineWidth', 1)
    scatter(1800, distdeg+0.15, 40, [0 1 0], 'filled', 'v', 'MarkerEdgeColor', 'k');
    text(1850, distdeg+0.1, 'NERSC2', 'Color', 'k', 'FontSize', 9)

    % NERSC3
    info_mooring = mooringinfo('NERSC3');
    stlo = info_mooring.rreflon;
    stla = info_mooring.rreflat;
    [~, wh] = min(abs(info_mooring.hmdep-1000));
    [t, x, jdn, num, depth] = getdvlaseis('NERSC3', info_mooring.snx(wh), ...
        '2020-01-08 12:07:00', '2020-01-08 12:37:00', fs);
    dt = datetime(jdn, 'ConvertFrom', 'datenum', ...
        'Format', 'uuuu-MM-dd''T''HH:mm:ss.SSSSSS', ...
        'TimeZone', 'UTC') + seconds(t);
    t_nersc1 = seconds(dt - dt_origin);
    x = bandpass(x, fs, 0.5, 5, 2, 1, 'butter', 'linear');
    [distkm, distdeg] = grcdist([evlo evla], [stlo stla]);
    plot(t_nersc1, x / max(abs(x)) + distdeg, ...
        'Color', 'k', 'LineWidth', 1)
    scatter(1800, distdeg+0.15, 40, [1 1 0], 'filled', 'v', 'MarkerEdgeColor', 'k');
    text(1850, distdeg+0.1, 'NERSC3', 'Color', 'k', 'FontSize', 9)

    % SIO1
    info_mooring = mooringinfo('SIO1');
    stlo = info_mooring.rreflon;
    stla = info_mooring.rreflat;
    [~, wh] = min(abs(info_mooring.hmdep-1000));
    [t, x, jdn, num, depth] = getdvlaseis('SIO1', info_mooring.snx(wh), ...
        '2020-01-08 12:15:00', '2020-01-08 12:37:00', fs);
    dt = datetime(jdn, 'ConvertFrom', 'datenum', ...
        'Format', 'uuuu-MM-dd''T''HH:mm:ss.SSSSSS', ...
        'TimeZone', 'UTC') + seconds(t);
    t_nersc1 = seconds(dt - dt_origin);
    x = bandpass(x, fs, 0.5, 5, 2, 1, 'butter', 'linear');
    [distkm, distdeg] = grcdist([evlo evla], [stlo stla]);
    plot(t_nersc1, x / max(abs(x)) + distdeg, ...
        'Color', 'k', 'LineWidth', 1)
    scatter(1800, distdeg+0.15, 40, [1 0.6 0.2], 'filled', 'v', 'MarkerEdgeColor', 'k');
    text(1850, distdeg+0.1, 'SIO1', 'Color', 'k', 'FontSize', 9)

    % SIO2
    info_mooring = mooringinfo('SIO2');
    stlo = info_mooring.rreflon;
    stla = info_mooring.rreflat;
    [~, wh] = min(abs(info_mooring.hmdep-1000));
    [t, x, jdn, num, depth] = getdvlaseis('SIO2', info_mooring.snx(wh), ...
        '2020-01-08 12:07:00', '2020-01-08 12:37:00', fs);
    dt = datetime(jdn, 'ConvertFrom', 'datenum', ...
        'Format', 'uuuu-MM-dd''T''HH:mm:ss.SSSSSS', ...
        'TimeZone', 'UTC') + seconds(t);
    t_nersc1 = seconds(dt - dt_origin);
    x = bandpass(x, fs, 0.5, 5, 2, 1, 'butter', 'linear');
    [distkm, distdeg] = grcdist([evlo evla], [stlo stla]);
    plot(t_nersc1, x / max(abs(x)) + distdeg, ...
        'Color', 'k', 'LineWidth', 1)
    scatter(1800, distdeg+0.15, 40, [0 1 1], 'filled', 'v', 'MarkerEdgeColor', 'k');
    text(1850, distdeg+0.1, 'SIO2', 'Color', 'k', 'FontSize', 9)

    % SIO3
    info_mooring = mooringinfo('SIO3');
    stlo = info_mooring.rreflon;
    stla = info_mooring.rreflat;
    [~, wh] = min(abs(info_mooring.hmdep-1000));
    [t, x, jdn, num, depth] = getdvlaseis('SIO3', info_mooring.snx(wh), ...
        '2020-01-08 12:07:00', '2020-01-08 12:37:00', fs);
    dt = datetime(jdn, 'ConvertFrom', 'datenum', ...
        'Format', 'uuuu-MM-dd''T''HH:mm:ss.SSSSSS', ...
        'TimeZone', 'UTC') + seconds(t);
    t_nersc1 = seconds(dt - dt_origin);
    x = bandpass(x, fs, 0.5, 5, 2, 1, 'butter', 'linear');
    [distkm, distdeg] = grcdist([evlo evla], [stlo stla]);
    plot(t_nersc1, x / max(abs(x)) + distdeg, ...
        'Color', 'k', 'LineWidth', 1)
    scatter(1800, distdeg+0.15, 40, [1 0.5 1], 'filled', 'v', 'MarkerEdgeColor', 'k');
    text(1850, distdeg+0.1, 'SIO3', 'Color', 'k', 'FontSize', 9)

    % add expected arrival time
    gcarcdists = (0:0.1:25)';
    [~, tP, tT1400] = PTtraveltime(gcarcdists, evdp, 0, 1.4);
    [~, ~, tT1600] = PTtraveltime(gcarcdists, evdp, 0, 1.6);
    plot(tP, gcarcdists, 'LineWidth', 1, 'Color', 'b')
    plot(tT1600, gcarcdists, 'LineWidth', 1, 'Color', [0 0.5 0], 'LineStyle', '--')
    plot(tT1400, gcarcdists, 'LineWidth', 1, 'Color', [0 0.5 0], 'LineStyle', '--')
    text(300, 18, 'P', 'FontSize', 10, 'Color', 'b')
    text(1320, 18, 'T', 'FontSize', 10, 'Color', [0 0.5 0])
    text(1200, 18, '1.6 km/s', 'FontSize', 10, 'Color', [0 0.5 0], 'Rotation', 60)
    text(1480, 18, '1.4 km/s', 'FontSize', 10, 'Color', [0 0.5 0], 'Rotation', 58)

    grid on
    xlim([0 2400])
    xticks(0:300:2400)
    xticklabels(0:5:40)
    xlabel('time since the origin (minutes)')
    ylabel('epicentral distance (degrees)')
    tit = title('Earthquake: M4.8 2020-01-08 12:07:46.392000, ID: 11167897');
    tit.Position(2) = 25.5;
    set(gca, 'Box', 'on', 'FontSize', 12, 'TickDir', 'out')
    set(gcf, 'Renderer', 'painters')
    figdisp(sprintf('%s_373_11167897', mfilename), [], [], 2, [], 'epstopdf')
elseif ev == 291
    fs = 50;

    ev291 = irisFetch.Events('MinimumMagnitude', 3, 'latitude', 69.575, ...
        'longitude', -144.472, 'maximumRadius', 1, ...
        'startTime', '2019-10-19T00:00:00', ...
        'endTime', '2019-10-19T00:59:59')';
    dt_origin = datetime(ev291.PreferredTime, ...
        'Format', 'dd-MMM-uuuu HH:mm:ss.SSSSSS', 'TimeZone', 'UTC');
    evlo = ev291.PreferredLongitude;
    evla = ev291.PreferredLatitude;
    evdp = ev291.PreferredDepth;
    mag = ev291.PreferredMagnitudeValue;
    
    % seismogram plot
    fig3 = figure(3);
    clf
    set(fig3, 'Units', 'inches', 'Position', [0 1 6 8])
    subplot('Position', [0.12 0.08 0.84 0.82])
    hold on

    % Land stations
    ddir = '/Users/spipatprathanporn/research/IFILES/SEISMOQUERY/Event_11133722/';

    station_list = {'II.ALE', 'TA.C24K', 'TA.D28M', 'IU.KBS'};
    
    for ii = 1:length(station_list)
        station = station_list{ii};
        search_str1 = sprintf('%s*BHE*.sac', station);
        search_str2 = sprintf('%s*BH2*.sac', station);
        try
            fname = cindeks(ls2cell(fullfile(ddir, search_str1), 1), 1);
        catch
            fname = cindeks(ls2cell(fullfile(ddir, search_str2), 1), 1);
        end
        [x, HdrData] = readsac(fname);
        [~, ~, ~, fs, ~, dt, ~] = gethdrinfo(HdrData);
        stlo = HdrData.STLO;
        stla = HdrData.STLA;
    
        sacpzfile = [fname, 'pz'];
        x = detrend(x, 2) .* shanning(length(x), 0.016);
        x = real(transfer(x, 1/fs, [0.01 0.02 8 10], resp_output_fmt, sacpzfile, 'sacpz'));
        x = bandpass(x, fs, 2, 6, 2, 1, 'butter', 'linear');
        t = seconds(dt - dt_origin);
        [distkm, distdeg] = grcdist([evlo evla], [stlo stla]);
        if distdeg < 20
            plot(t, x / max(abs(x)) * 2 + distdeg, 'Color', [0.65 0.65 0.65], 'LineWidth', 1)
        else
            plot(t, x / max(abs(x)) + distdeg, 'Color', [0.65 0.65 0.65], 'LineWidth', 1)
        end
        text(1800, distdeg+0.40, station, 'FontSize', 9, ...
            'Color', [0.3 0.3 0.3])
    end

    % % NERSC1
    % info_mooring = mooringinfo('NERSC1');
    % stlo = info_mooring.rreflon;
    % stla = info_mooring.rreflat;
    % [~, wh] = min(abs(info_mooring.hmdep-1000));
    % [t, x, jdn, num, depth] = getdvlaseis('NERSC1', info_mooring.snx(wh), ...
    %     '2019-10-19 00:01:00', '2019-10-19 00:35:00', fs);
    % dt = datetime(jdn, 'ConvertFrom', 'datenum', ...
    %     'Format', 'uuuu-MM-dd''T''HH:mm:ss.SSSSSS', ...
    %     'TimeZone', 'UTC') + seconds(t);
    % t_nersc1 = seconds(dt - dt_origin);
    % x = bandpass(x, fs, 0.5, 5, 2, 1, 'butter', 'linear');
    % [distkm, distdeg] = grcdist([evlo evla], [stlo stla]);
    % plot(t_nersc1, x / max(abs(x)) + distdeg, ...
    %     'Color', 'k', 'LineWidth', 1)
    % scatter(2100, distdeg+0.15, 40, [1 0 0], 'filled', 'v', 'MarkerEdgeColor', 'k');
    % text(2150, distdeg+0.1, 'NERSC1', 'Color', 'k', 'FontSize', 9)
    % 
    % % NERSC2
    % info_mooring = mooringinfo('NERSC2');
    % stlo = info_mooring.rreflon;
    % stla = info_mooring.rreflat;
    % [~, wh] = min(abs(info_mooring.hmdep-1000));
    % [t, x, jdn, num, depth] = getdvlaseis('NERSC2', info_mooring.snx(wh), ...
    %     '2019-10-19 00:01:00', '2019-10-19 00:35:00', fs);
    % dt = datetime(jdn, 'ConvertFrom', 'datenum', ...
    %     'Format', 'uuuu-MM-dd''T''HH:mm:ss.SSSSSS', ...
    %     'TimeZone', 'UTC') + seconds(t);
    % t_nersc1 = seconds(dt - dt_origin);
    % x = bandpass(x, fs, 0.5, 5, 2, 1, 'butter', 'linear');
    % [distkm, distdeg] = grcdist([evlo evla], [stlo stla]);
    % plot(t_nersc1, x / max(abs(x)) + distdeg, ...
    %     'Color', 'k', 'LineWidth', 1)
    % scatter(2100, distdeg+0.15, 40, [0 1 0], 'filled', 'v', 'MarkerEdgeColor', 'k');
    % text(2150, distdeg+0.1, 'NERSC2', 'Color', 'k', 'FontSize', 9)

    % NERSC3
    info_mooring = mooringinfo('NERSC3');
    stlo = info_mooring.rreflon;
    stla = info_mooring.rreflat;
    [~, wh] = min(abs(info_mooring.hmdep-1000));
    [t, x, jdn, num, depth] = getdvlaseis('NERSC3', info_mooring.snx(end), ...
        '2019-10-19 00:01:00', '2019-10-19 00:45:00', fs);
    dt = datetime(jdn, 'ConvertFrom', 'datenum', ...
        'Format', 'uuuu-MM-dd''T''HH:mm:ss.SSSSSS', ...
        'TimeZone', 'UTC') + seconds(t);
    t_nersc1 = seconds(dt - dt_origin);
    x = bandpass(x, fs, 5, 8, 2, 1, 'butter', 'linear');
    [distkm, distdeg] = grcdist([evlo evla], [stlo stla]);
    plot(t_nersc1, x / max(abs(x)) + distdeg, ...
        'Color', 'k', 'LineWidth', 1)
    scatter(2100, distdeg+1.15, 40, [1 1 0], 'filled', 'v', 'MarkerEdgeColor', 'k');
    text(2150, distdeg+1.1, 'NERSC3', 'Color', 'k', 'FontSize', 9)

    % SIO1
    info_mooring = mooringinfo('SIO1');
    stlo = info_mooring.rreflon;
    stla = info_mooring.rreflat;
    [~, wh] = min(abs(info_mooring.hmdep-1000));
    [t, x, jdn, num, depth] = getdvlaseis('SIO1', info_mooring.snx(wh), ...
        '2019-10-19 00:01:00', '2019-10-19 00:35:00', fs);
    dt = datetime(jdn, 'ConvertFrom', 'datenum', ...
        'Format', 'uuuu-MM-dd''T''HH:mm:ss.SSSSSS', ...
        'TimeZone', 'UTC') + seconds(t);
    t_nersc1 = seconds(dt - dt_origin);
    x = bandpass(x, fs, 5, 8, 2, 1, 'butter', 'linear');
    [distkm, distdeg] = grcdist([evlo evla], [stlo stla]);
    plot(t_nersc1, x / max(abs(x)) + distdeg, ...
        'Color', 'k', 'LineWidth', 1)
    scatter(2100, distdeg+0.15, 40, [1 0.6 0.2], 'filled', 'v', 'MarkerEdgeColor', 'k');
    text(2150, distdeg+0.1, 'SIO1', 'Color', 'k', 'FontSize', 9)

    % SIO2
    info_mooring = mooringinfo('SIO2');
    stlo = info_mooring.rreflon;
    stla = info_mooring.rreflat;
    [~, wh] = min(abs(info_mooring.hmdep-1000));
    [t, x, jdn, num, depth] = getdvlaseis('SIO2', info_mooring.snx(wh), ...
        '2019-10-19 00:01:00', '2019-10-19 00:35:00', fs);
    dt = datetime(jdn, 'ConvertFrom', 'datenum', ...
        'Format', 'uuuu-MM-dd''T''HH:mm:ss.SSSSSS', ...
        'TimeZone', 'UTC') + seconds(t);
    t_nersc1 = seconds(dt - dt_origin);
    x = bandpass(x, fs, 5, 8, 2, 1, 'butter', 'linear');
    [distkm, distdeg] = grcdist([evlo evla], [stlo stla]);
    plot(t_nersc1, x / max(abs(x)) + distdeg, ...
        'Color', 'k', 'LineWidth', 1)
    scatter(2100, distdeg+0.15, 40, [0 1 1], 'filled', 'v', 'MarkerEdgeColor', 'k');
    text(2150, distdeg+0.1, 'SIO2', 'Color', 'k', 'FontSize', 9)

    % SIO3
    info_mooring = mooringinfo('SIO3');
    stlo = info_mooring.rreflon;
    stla = info_mooring.rreflat;
    [~, wh] = min(abs(info_mooring.hmdep-1000));
    [t, x, jdn, num, depth] = getdvlaseis('SIO3', info_mooring.snx(wh), ...
        '2019-10-19 00:01:00', '2019-10-19 00:35:00', fs);
    dt = datetime(jdn, 'ConvertFrom', 'datenum', ...
        'Format', 'uuuu-MM-dd''T''HH:mm:ss.SSSSSS', ...
        'TimeZone', 'UTC') + seconds(t);
    t_nersc1 = seconds(dt - dt_origin);
    x = bandpass(x, fs, 5, 8, 2, 1, 'butter', 'linear');
    [distkm, distdeg] = grcdist([evlo evla], [stlo stla]);
    plot(t_nersc1, x / max(abs(x)) + distdeg, ...
        'Color', 'k', 'LineWidth', 1)
    scatter(2100, distdeg+0.15, 40, [1 0.5 1], 'filled', 'v', 'MarkerEdgeColor', 'k');
    text(2150, distdeg+0.1, 'SIO3', 'Color', 'k', 'FontSize', 9)

    % add expected arrival time
    gcarcdists = (0:0.1:35)';
    [~, tP, tT1400] = PTtraveltime(gcarcdists, evdp, 0, 1.4);
    [~, ~, tT1600] = PTtraveltime(gcarcdists, evdp, 0, 1.6);
    plot(tP, gcarcdists, 'LineWidth', 1, 'Color', 'b')
    plot(tT1600, gcarcdists, 'LineWidth', 1, 'Color', [0 0.5 0], 'LineStyle', '--')
    plot(tT1400, gcarcdists, 'LineWidth', 1, 'Color', [0 0.5 0], 'LineStyle', '--')
    text(300, 18, 'P', 'FontSize', 10, 'Color', 'b')
    text(1320, 18, 'T', 'FontSize', 10, 'Color', [0 0.5 0])
    text(1200, 18, '1.6 km/s', 'FontSize', 10, 'Color', [0 0.5 0], 'Rotation', 54.5)
    text(1480, 18, '1.4 km/s', 'FontSize', 10, 'Color', [0 0.5 0], 'Rotation', 52.5)

    grid on
    xlim([0 2400])
    ylim([0 35])
    xticks(0:300:2400)
    xticklabels(0:5:40)
    xlabel('time since the origin (minutes)')
    ylabel('epicentral distance (degrees)')
    tit = title('Earthquake: M3.3 2019-10-19 00:01:13.158000, ID: 11133722');
    tit.Position(2) = 35.5;
    set(gca, 'Box', 'on', 'FontSize', 12, 'TickDir', 'out')
    set(gcf, 'Renderer', 'painters')
    figdisp(sprintf('%s_291_11133722', mfilename), [], [], 2, [], 'epstopdf')

% the uncatalogged event
elseif ev == 294
    ax = polarbasemap(15, 'notcartoon');
    caatex = polaraddcaatex(ax);

    % longitude label
    for ii = 0:30:330
        if and(ii>0, ii<180)
            plotpolarmap(@text, ii, 74.5, ...
                [sprintf('%d', ii) '^{\circ}'], 'FontSize', 12, ...
                'HorizontalAlignment', 'left');
        elseif and(ii>180, ii<360)
            plotpolarmap(@text, ii, 74.5, ...
                [sprintf('%d', ii) '^{\circ}'], 'FontSize', 12, ...
                'HorizontalAlignment', 'right');
        else
            plotpolarmap(@text, ii, 74.5, ...
                [sprintf('%d', ii) '^{\circ}'], 'FontSize', 12, ...
                'HorizontalAlignment', 'center');
        end
    end

    tP = [1855, 1865, 1857];
    tT = [2000, 2050, 2010];
    

    tT_SIO1 = 2015.5;

    dists = zeros(size(tP));
    for ii = 1:3
        dists(ii) =  PTtdiff2distance(tT(ii)-tP(ii), 10, 0, 1.5);
    end

    % NERSC1-3
    stla = [84.0033, 82.5125, 83.4420]';
    stlo = [28.3916, 23.9628, 25.6672]';
    stco = [1 0 0; 0 0.75 0; 1 1 0];

    % draw the circles
    theta = (0:10:360)';
    l = gobjects(3,1);
    for ii = 1:3
        lats = nan(size(theta));
        lons = nan(size(theta));
        for jj = 1:length(theta)
            [lats(jj), lons(jj)] = reckon(stla(ii), stlo(ii), ...
                dists(ii), theta(jj));
        end
        l(ii) = plotpolarmap(@plot, lons, lats, 'LineWidth', 2, ...
                'Color', stco(ii,:));
    end

    % draw the SIO1 circle
    dist_SIO1 = 14.12;
    stlo_SIO1 = -148.3434;
    stla_SIO1 = 80.8697;
    for jj = 1:length(theta)
        [lats(jj), lons(jj)] = reckon(stla_SIO1, stlo_SIO1, ...
            dist_SIO1, theta(jj));
        
    end
    plotpolarmap(@plot, lons, lats, 'LineWidth', 2, 'Color', [1 0.5 0]);
end


end