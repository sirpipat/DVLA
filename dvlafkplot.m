function dvlafkplot(station, datadir, metadatadir)
% DVLAFKPLOT(station, masterdir, metadatadir)
%
% Make F-K plots from an entire array at a station from CAATEX Experiment
%
% INPUT:
% station       station name: either 'SIO1', 'SIO2', 'SIO3', 'NERSC1',
%               'NERSC2', or 'NERSC3'
% datadir       directory containing data from all stations
% metadatadir   directory containing metadata for the station
%
% Last modified by spipatprathanporn@ucsd.edu, 04/06/2026

% TODO: change to environment variables
% TODO: decide where to host data+metadata
defval('metadatadir', '/Volumes/aog/AOGquake/work/Ivy/CAATEX/metadata');
defval('datadir', '/Volumes/aog/AOGquake/vault/CAATEX/')
addpath(genpath(metadatadir))

% determine which directories / metadata files to read
if strcmpi(station, 'SIO1')
    masterdir = fullfile(datadir, 'SIO1', 'HMrcv/');
    metadatafile = fullfile(metadatadir, 'minfo', 'minfo.SN101.mat');
    load(metadatafile, 'minfo');
    snx = minfo.snx;
    hmdep = minfo.hmdep;
    metadatafile = fullfile(metadatadir, 'minfo', 'minfo.SN102.mat');
    load(metadatafile, 'minfo');
    snx = [snx minfo.snx];
    hmdep = [hmdep minfo.hmdep];
elseif strcmpi(station, 'SIO2')
    masterdir = fullfile(datadir, 'SIO2', 'HMrcv/');
    metadatafile = fullfile(metadatadir, 'minfo', 'minfo.SN103.mat');
    load(metadatafile, 'minfo');
    snx = minfo.snx;
    hmdep = minfo.hmdep;
    metadatafile = fullfile(metadatadir, 'minfo', 'minfo.SN104.mat');
    load(metadatafile, 'minfo');
    snx = [snx minfo.snx];
    hmdep = [hmdep minfo.hmdep];
elseif strcmpi(station, 'SIO3')
    masterdir = fullfile(datadir, 'SIO3', 'HMrcv/');
    metadatafile = fullfile(metadatadir, 'minfo', 'minfo.SN105.mat');
    load(metadatafile, 'minfo');
    snx = minfo.snx;
    hmdep = minfo.hmdep;
    metadatafile = fullfile(metadatadir, 'minfo', 'minfo.SN106.mat');
    load(metadatafile, 'minfo');
    snx = [snx minfo.snx];
    hmdep = [hmdep minfo.hmdep];
elseif strcmpi(station, 'NERSC1')
    masterdir = fullfile(datadir, 'NERSC1', 'HMrcv/');
    metadatafile = fullfile(metadatadir, 'minfo', 'minfo.SN110.mat');
    load(metadatafile, 'minfo');
    snx = minfo.snx;
    hmdep = minfo.hmdep;
elseif strcmpi(station, 'NERSC2')
    masterdir = fullfile(datadir, 'NERSC2', 'HMrcv/');
    metadatafile = fullfile(metadatadir, 'minfo', 'minfo.SN131.mat');
    load(metadatafile, 'minfo');
    snx = minfo.snx;
    hmdep = minfo.hmdep;
elseif strcmpi(station, 'NERSC3')
    masterdir = fullfile(datadir, 'NERSC3', 'HMrcv/');
    metadatafile = fullfile(metadatadir, 'minfo', 'minfo.SN133.mat');
    load(metadatafile, 'minfo');
    snx = minfo.snx;
    hmdep = minfo.hmdep;
else
    error('Invalid station name')
    return
end

% make F-K plot when first hydrophone can record acoustic signal
[allfiles, fndex] = allfile(sprintf('%s/%03d/', masterdir, snx(1)));

% loop over all sections
for ii = 1:fndex
    fname = removepath(allfiles{ii});
    
    % construct beginning and ending datetimes
    jday = str2double(fname(5:7));
    jhour = str2double(fname(8:9));
    jminute = str2double(fname(10:11));
    jsecond = str2double(fname(12:13));
    dt_begin = datetime(2018,12,31,jhour,jminute,jsecond+1, ...
        'Format', 'uuuu-MM-dd''T''HH:mm:ss.SSSSSS') + days(jday);
    dt_end = dt_begin + seconds(2710-1);

    % request data from all hydrophones from the mooring
    [t, xx, jdn, num, depth] = getdvlaseis(station, 'all', dt_begin, ...
        dt_end, 100, datadir, metadatadir);
    dt = t(2) - t(1);
    dz = hmdep(2) - hmdep(1);
    xx = xx';

    % F-K transformation
    [XX, F, K] = fk(xx, dt, dz);

    % plot
    figure(1);
    set(gcf, 'Unit', 'inches', 'Position', [0 1 11 5]);
    clf;
    ax1 = subplot('Position', [0.06 0.12 0.41 0.8]);
    imagesc(K, F, 10*log10(abs(XX).^2));
    grid on
    axis xy
    xlabel('wavenumber (m^{-1})')
    ylabel('frequency (Hz)')
    title(sprintf('station: %s - %s', station, dt_begin))
    cb = colorbar;
    cb.TickDirection = 'out';
    cb.Label.String = '10 log_{10}(spectral density) (dB)';
    set(gca, 'TickDir', 'out', 'Box', 'on', 'FontSize', 12)

    ax2 = subplot('Position', [0.56 0.12 0.41 0.8]);
    imagesc(K, F, 10*log10(abs(XX).^2));
    grid on
    axis xy
    xlabel('wavenumber (m^{-1})')
    ylabel('frequency (Hz)')
    title('zoomed-in')
    ylim([-5 5])
    cb = colorbar;
    cb.TickDirection = 'out';
    cb.Label.String = '10 log_{10}(spectral density) (dB)';
    set(gca, 'TickDir', 'out', 'Box', 'on', 'FontSize', 12)

    set(gcf, 'Renderer', 'painters')
    savename = sprintf('%s_%s_%s', mfilename, station, ...
        replace(fname, 'nc', 'eps'));
    figdisp(savename, [], [], 2, [], 'epstopdf')

    % make seismogram plot
    figure(2)
    set(gcf, 'Unit', 'inches', 'Position', [0 6 11 5]);
    clf
    ax21 = subplot('Position', [0.06 0.12 0.41 0.8]);
    hold on
    scalingfactor = indeks(maxk(max(abs(xx(:,9:end)), [], 1), 2), 2) / dz * 2;
    for kk = 1:size(xx, 2)
        plot(t, xx(:,kk) / scalingfactor + hmdep(kk), 'LineWidth', 1, 'Color', 'k')
    end
    xlabel('time (s)')
    ylabel('depth (m)')
    title(sprintf('station: %s - %s', station, dt_begin))
    grid on
    axis ij
    axis tight
    ylim([hmdep(1)-dz hmdep(end)+dz])
    set(ax21, 'FontSize', 12, 'TickDir', 'out', 'Box', 'on')

    ax22 = subplot('Position', [0.56 0.12 0.41 0.8]);
    hold on
    xxf = xx;
    for kk = 1:size(xx, 2)
        xxf(:,kk) = bandpass(detrend(xx(:,kk)) .* ...
            shanning(size(xx,1), 0.01), 1/dt, 0.5, 5, 2, 1, 'butter', ...
            'linear');
    end
    scalingfactor = indeks(maxk(max(abs(xxf(:,9:end)), [], 1), 2), 2) / dz * 1.5;
    for kk = 1:size(xx, 2)
        plot(t, xxf(:,kk) / scalingfactor + hmdep(kk), 'LineWidth', 1, 'Color', 'k')
    end
    xlabel('time (s)')
    ylabel('depth (m)')
    title('0.5--5 Hz')
    grid on
    axis ij
    axis tight
    ylim([hmdep(1)-dz hmdep(end)+dz])
    set(ax22, 'FontSize', 12, 'TickDir', 'out', 'Box', 'on')

    set(gcf, 'Renderer', 'painters')
    savename = sprintf('%s_%s_%s_%s', mfilename, 't-series', station, ...
        replace(fname, 'nc', 'eps'));
    figdisp(savename, [], [], 2, [], 'epstopdf')
end
end