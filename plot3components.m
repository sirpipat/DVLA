function [st, hdr] = plot3components(ddir, net, sta, loc, ev, plt)
% [st, hdr] = PLOT3COMPONENTS(ddir, net, sta, loc, ev, plt)
%
% Plot 3-component seismograms at a given station. Then, mark the expected
% arrival times of P, PP, PcP, S, SS, ScS, T-phase arrival times given the
% earthquake event.
%
% INPUT:
% ddir          directory to the SAC files
% net           network code. Use '*' for wildcards
% sta           station code. Use '*' for wildcards
% loc           location. Use '*' for wildcards
% ev            an event struct with a following fields
%       - PreferredTime             'uuuu-MM-dd hh:mm:ss.SSSSSS'
%       - PreferredLatitude
%       - PreferredLongitude
%       - PreferredDepth            in km
%       - PreferredMagnitudeValue
%       - PreferredMagnitudeType
%       - PublicId                  'query?eventid=xxxxxxxx'
% plt           whether to plot or not
%
% OUTPUT:
% st            three component seismograms
% hdr           seismograms' metadata headers
%
% Last modified by spipatprathanporn@ucsd.edu, 03/25/2026

defval('net', '*')
defval('sta', '*')
defval('loc', '*')
defval('plt', true)

% validate input
net = strip(net);
sta = strip(sta);
loc = strip(loc);

dt0 = datetime(ev.PreferredTime, 'TimeZone', 'UTC');

resp_output_fmt = 'velocity';
passband = [2 8];

sacfilesZ = ls2cell(fullfile(ddir, sprintf('%s.%s.%s.BHZ.*.sac', ...
    net, sta, loc)), 1);
if ~isempty(sacfilesZ)
    sacfileZ = sacfilesZ{1};
    [SeisZ, HdrData] = readsac(sacfileZ);
    [~, ~, ~, fs, ~, dts, ~] = gethdrinfo(HdrData);
    tZ = seconds(dts - dt0);
    try
        sacpzfileZ = [sacfileZ 'pz'];
        SeisZ = detrend(SeisZ, 2) .* shanning(length(SeisZ), 0.016);
        SeisZ = real(transfer(SeisZ, 1/fs, [0.01 0.02 8 10], ...
            resp_output_fmt, sacpzfileZ, 'sacpz'));
        st{1} = [tZ, SeisZ];
        HdrData.KCMPNM = 'BHZ';
        hdr{1} = update_header_event_info(HdrData, ev);
    catch ME
        ME.getReport()
    end

    sacfileN = replace(sacfileZ, 'BHZ', 'BHN');
    try
        [SeisN, HdrData] = readsac(sacfileN);
        [~, ~, ~, fs, ~, dts, ~] = gethdrinfo(HdrData);
        tN = seconds(dts - dt0);
        sacpzfileN = [sacfileN, 'pz'];
        SeisN = detrend(SeisN, 2) .* shanning(length(SeisN), 0.016);
        SeisN = real(transfer(SeisN, 1/fs, [0.01 0.02 8 10], ...
            resp_output_fmt, sacpzfileN, 'sacpz'));
        st{2} = [tN, SeisN];
        HdrData.KCMPNM = 'BHN';
        hdr{2} = update_header_event_info(HdrData, ev);

        sacfileE = replace(sacfileZ, 'BHZ', 'BHE');
        [SeisE, HdrData] = readsac(sacfileE);
        [~, ~, ~, fs, ~, dts, ~] = gethdrinfo(HdrData);
        tE = seconds(dts - dt0);
        sacpzfileE = [sacfileE, 'pz'];
        SeisE = detrend(SeisE, 2) .* shanning(length(SeisE), 0.016);
        SeisE = real(transfer(SeisE, 1/fs, [0.01 0.02 8 10], ...
            resp_output_fmt, sacpzfileE, 'sacpz'));
        st{3} = [tE, SeisE];
        HdrData.KCMPNM = 'BHE';
        hdr{3} = update_header_event_info(HdrData, ev);
    catch ME
        % Either N or E component file do not exist
        % Search for 1,2 component instead
        if contains(ME.message, 'does not exist')
            try
                sacfile1 = replace(sacfileZ, 'BHZ', 'BH1');
                [Seis1, HdrData] = readsac(sacfile1);
                [~, ~, ~, fs, ~, dts, ~] = gethdrinfo(HdrData);
                t1 = seconds(dts - dt0);
                sacpzfile1 = [sacfile1, 'pz'];
                Seis1 = detrend(Seis1, 2) .* shanning(length(Seis1), 0.016);
                Seis1 = real(transfer(Seis1, 1/fs, [0.01 0.02 8 10], ...
                    resp_output_fmt, sacpzfile1, 'sacpz'));
                st{2} = [t1, Seis1];
                HdrData.KCMPNM = 'BH1';
                hdr{2} = update_header_event_info(HdrData, ev);

                sacfile2 = replace(sacfileZ, 'BHZ', 'BH2');
                [Seis2, HdrData] = readsac(sacfile2);
                [~, ~, ~, fs, ~, dts, ~] = gethdrinfo(HdrData);
                t2 = seconds(dts - dt0);
                sacpzfile2 = [sacfile2, 'pz'];
                Seis2 = detrend(Seis2, 2) .* shanning(length(Seis2), 0.016);
                Seis2 = real(transfer(Seis2, 1/fs, [0.01 0.02 8 10], ...
                    resp_output_fmt, sacpzfile2, 'sacpz'));
                st{3} = [t2, Seis2];
                HdrData.KCMPNM = 'BH2';
                hdr{3} = update_header_event_info(HdrData, ev);
            catch ME
                ME.getReport()
            end
        else
            ME.getReport()
        end
    end
else
    try
        sacfile3 = ls2cell(fullfile(ddir, sprintf('%s.%s.%s.BH3.*.%s.sac', ...
            ddir, net, sta, loc)), 1);
        [Seis3, HdrData] = readsac(sacfile3);
        [~, ~, ~, fs, ~, dts, ~] = gethdrinfo(HdrData);
        t3 = seconds(dts - dt0);
        sacpzfile3 = [sacfile3, 'pz'];
        Seis3 = detrend(Seis3, 2) .* shanning(length(Seis3), 0.016);
        Seis3 = real(transfer(Seis3, 1/fs, [0.01 0.02 8 10], ...
            resp_output_fmt, sacpzfile3, 'sacpz'));
        st{1} = [t3, Seis3];
        HdrData.KCMPNM = 'BH3';
        hdr{1} = update_header_event_info(HdrData, ev);

        sacfile1 = replace(sacfile3, 'BH3', 'BH1');
        [Seis1, HdrData] = readsac(sacfile1);
        [~, ~, ~, fs, ~, dts, ~] = gethdrinfo(HdrData);
        t1 = seconds(dts - dt0);
        sacpzfile1 = [sacfile1, 'pz'];
        Seis1 = detrend(Seis1, 2) .* shanning(length(Seis1), 0.016);
        Seis1 = real(transfer(Seis1, 1/fs, [0.01 0.02 8 10], ...
            resp_output_fmt, sacpzfile1, 'sacpz'));
        st{2} = [t1, Seis1];
        HdrData.KCMPNM = 'BH1';
        hdr{2} = update_header_event_info(HdrData, ev);

        sacfile2 = replace(sacfile3, 'BH3', 'BH2');
        [Seis2, HdrData] = readsac(sacfile2);
        [~, ~, ~, fs, ~, dts, ~] = gethdrinfo(HdrData);
        t2 = seconds(dts - dt0);
        sacpzfile2 = [sacfile2, 'pz'];
        Seis2 = detrend(Seis2, 2) .* shanning(length(Seis2), 0.016);
        Seis2 = real(transfer(Seis2, 1/fs, [0.01 0.02 8 10], ...
            resp_output_fmt, sacpzfile2, 'sacpz'));
        st{3} = [t2, Seis2];
        HdrData.KCMPNM = 'BH2';
        hdr{3} = update_header_event_info(HdrData, ev);

        catch ME
            ME.getReport()
    end
end

% plot
if plt
    % compute travel times
    tP = getfield(indeks(tauptime('dep', hdr{1}.EVDP, 'ph', 'p,P', ...
        'deg', hdr{1}.GCARC, 'mod', 'ak135'), 1), 'time');
    tS = getfield(indeks(tauptime('dep', hdr{1}.EVDP, 'ph', 's,S', ...
        'deg', hdr{1}.GCARC, 'mod', 'ak135'), 1), 'time');
    tPP = getfield(indeks(tauptime('dep', hdr{1}.EVDP, 'ph', 'PP', ...
        'deg', hdr{1}.GCARC, 'mod', 'ak135'), 1), 'time');
    tSS = getfield(indeks(tauptime('dep', hdr{1}.EVDP, 'ph', 'SS', ...
        'deg', hdr{1}.GCARC, 'mod', 'ak135'), 1), 'time');
    tPcP = getfield(indeks(tauptime('dep', hdr{1}.EVDP, 'ph', 'PcP', ...
        'deg', hdr{1}.GCARC, 'mod', 'ak135'), 1), 'time');
    tScS = getfield(indeks(tauptime('dep', hdr{1}.EVDP, 'ph', 'ScS', ...
        'deg', hdr{1}.GCARC, 'mod', 'ak135'), 1), 'time');
    [~,~,tT1600] = PTtraveltime(hdr{1}.GCARC, hdr{1}.EVDP, 0, 1.6, 0);
    [~,~,tT1400] = PTtraveltime(hdr{1}.GCARC, hdr{1}.EVDP, 0, 1.4, 0);

    % compute phase arrivals
    fig = figure(10);
    clf
    set(gcf, 'Units', 'inches', 'Position', [0 1 10 5])

    for ii = 1:length(st)
        subplot(3,1,ii)
        xf = bandpass(st{ii}(:,2), 1/hdr{ii}.DELTA, passband(1), ...
            passband(2), 2, 1, 'butter', 'linear');
        plot(st{ii}(:,1), xf, 'LineWidth', 1, 'Color', 'k')
        grid on
        axis tight
        hold on
        xlim([tP-50 tT1400+500])
        ylim([-1.15 1.15] * max(abs(xf(and(st{ii}(:,1) > tP-50, ...
            st{ii}(:,1) < tT1400+500)))))
        vline(gca, tP, 'LineWidth', 1, 'Color', [0 0 1]);
        vline(gca, tPP, 'LineWidth', 1, 'Color', [0.1 0.3 1]);
        vline(gca, tPcP, 'LineWidth', 1, 'Color', [0.4 0.6 1]);
        vline(gca, tS, 'LineWidth', 1, 'Color', [1 0 0]);
        vline(gca, tSS, 'LineWidth', 1, 'Color', [1 0.3 0.1]);
        vline(gca, tScS, 'LineWidth', 1, 'Color', [1 0.6 0.4]);
        vline(gca, tT1600, 'LineWidth', 1, 'Color', [0 0.6 0]);
        vline(gca, tT1400, 'LineWidth', 1, 'Color', [0 0.6 0]);

        set(gca, 'Box', 'on', 'TickDir', 'out', 'FontSize', 12)
        
        title(sprintf('%s.%s (%.2f deg, %d km): %s', strip(hdr{ii}.KNETWK), ...
            strip(hdr{ii}.KSTNM), hdr{ii}.GCARC, round(deg2km(hdr{ii}.GCARC)), ...
            strip(hdr{ii}.KCMPNM)));

        if ii == length(st)
            xlabel('time since origin (s)')
        else
            nolabels(gca, 1)
        end
    end

    set(gcf, 'Renderer', 'painters')

    savename = sprintf('%s_%d_%s.%s.%s.eps', mfilename, hdr{1}.USER7, ...
        strip(hdr{1}.KNETWK), strip(hdr{1}.KSTNM), strip(hdr{1}.KHOLE));
    figdisp(savename, [], [], 2, [], 'epstopdf')
    delete(fig)
end
end

function hdr = update_header_event_info(hdr, ev)
% include event info to the SAC header
dt0 = datetime(ev.PreferredTime, 'TimeZone', 'UTC');
evlo = ev.PreferredLongitude;
evla = ev.PreferredLatitude;
evdp = ev.PreferredDepth;
mag = ev.PreferredMagnitudeValue;
eventid = str2double(indeks(split(ev.PublicId, '='), 'end'));

hdr.EVLO = evlo;
hdr.EVLA = evla;
hdr.EVDP = evdp;
hdr.MAG = mag;
hdr.USER7 = eventid;

% compute the great-circle epicentral distance
[azdeg, bazdeg, distdeg] = azimdist([evlo evla], [hdr.STLO hdr.STLA]);
hdr.AZ = azdeg;
hdr.BAZ = bazdeg;
hdr.GCARC = distdeg;

% origin time relative to seismogram's reference time in seconds is
% recorded in a variable name USER8
dt_ref = gethdrinfo(hdr);
hdr.USER8 = seconds(dt0 -dt_ref);
end