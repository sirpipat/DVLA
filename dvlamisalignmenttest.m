function dvlamisalignmenttest(station, fs, metadatadir, datadir)
% DVLAMISALIGNMENTTEST(station, fs)
%
% Find any time misalignment between the top and bottom arrays at the
% hydrophone mooring. It plots the timeshift between two adjacent
% hydrophones in the same array (19-20) and two adjacent hydrophones in
% different arrays (20-21).
%
% INPUT:
% station       station name either 'SIO1', 'SIO2', or 'SIO3'
% fs            sampling rate
%
% Last modified by spipatprathanporn@ucsd.edu, 03/04/2026

% TODO: change to environment variables
% TODO: decide where to host data+metadata
defval('metadatadir', '/Volumes/aog/AOGquake/work/Ivy/CAATEX/metadata');
defval('datadir', '/Volumes/aog/AOGquake/vault/CAATEX/')
defval('fs', [])

% get metadata for the mooring
minfo = mooringinfo(station, metadatadir);
snx = minfo.snx;

% determine which directories to read
masterdir = fullfile(datadir, station, 'HMrcv');

% second to the bottom of the top array (for comparison)
[allfiles19, sndex19] = allfile(sprintf('%s/%03d/', masterdir, snx(19)));

% collecting day
d = nan([sndex19 1]);

% collecting timeshift parameters
t19_20 = nan([sndex19 1]);
t20_21 = nan([sndex19 1]);

for ii = 1:sndex19
    try
        [x19, t19, jdn19] = HMrcvx(allfiles19{ii});
        [x20, t20, jdn20] = HMrcvx(replace(allfiles19{ii}, ...
            sprintf('/%03d/', snx(19)), ...
            sprintf('/%03d/', snx(20))));
        [x21, t21, jdn21] = HMrcvx(replace(allfiles19{ii}, ...
            sprintf('/%03d/', snx(19)), ...
            sprintf('/%03d/', snx(21))));
    catch ME
        ME.getReport()
        continue
    end

    % lowpass filter
    fs19 = 1 / (t19(2) - t19(1));
    x19 = lowpass(x19, fs19, fs/2, 2, 1, 'butter', 'linear');
    fs20 = 1 / (t20(2) - t20(1));
    x20 = lowpass(x20, fs20, fs/2, 2, 1, 'butter', 'linear');
    fs21 = 1 / (t21(2) - t21(1));
    x21 = lowpass(x21, fs21, fs/2, 2, 1, 'butter', 'linear');

    % interpolate to the same time
    dt_ref19 = datetime(jdn19, 'ConvertFrom', 'datenum', ...
        'Format', 'uuuu-MM-dd''T''HH:mm:ss.SSSSSS');
    dt19 = dt_ref19 + seconds(t19);

    dt_ref20 = datetime(jdn20, 'ConvertFrom', 'datenum', ...
        'Format', 'uuuu-MM-dd''T''HH:mm:ss.SSSSSS');
    dt20 = dt_ref20 + seconds(t20);

    dt_ref21 = datetime(jdn21, 'ConvertFrom', 'datenum', ...
        'Format', 'uuuu-MM-dd''T''HH:mm:ss.SSSSSS');
    dt21 = dt_ref21 + seconds(t21);

    % time to convert
    dt = datetime(ceil(jdn20*2)/2, 'ConvertFrom', 'datenum', ...
        'Format', 'uuuu-MM-dd''T''HH:mm:ss.SSSSSS') + ...
        seconds(0:(1/fs):2700);

    x19int = interp1(dt19, x19, dt, 'linear');
    x20int = interp1(dt20, x20, dt, 'linear');
    x21int = interp1(dt21, x21, dt, 'linear');

    % cross-correlate (limit max lag to 10 s)
    [r19_20, lag19_20] = xcorr(x19int, x20int, 10*fs, 'coeff');
    [r20_21, lag20_21] = xcorr(x20int, x21int, 10*fs, 'coeff');


    t19_20(ii) = lag19_20(r19_20 == max(r19_20)) / fs;
    t20_21(ii) = lag20_21(r20_21 == max(r20_21)) / fs;
    d(ii) = days(dt(1) - datetime(2018, 12, 31, 0, 0, 0, 0, ...
        'Format', 'uuuu-MM-dd''T''HH:mm:ss.SSSSSS'));
end

figure(3)
clf
set(gcf, 'Units', 'inches', 'Position', [0 1 8 6]);
hold on
plot(d, t19_20, 'LineWidth', 1, 'Color', 'k')
plot(d, t20_21, 'LineWidth', 1, 'Color', [0.8 0.2 0.2])
grid on
xlabel('day')
ylabel('optimal timeshift (s)')
title(sprintf('%s: fs = %.2f Hz', station, fs))
legend('same array (19-20)', 'different arrays (20-21)', 'Location', 'southwest')
set(gca, 'FontSize', 12, 'TickDir', 'out', 'Box', 'on')
set(gcf, 'Renderer', 'painters')
end