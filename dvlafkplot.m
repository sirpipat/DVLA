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
% Last modified by spipatprathanporn@ucsd.edu, 12/02/2025

% TODO: change to environment variables
% TODO: decide where to host data+metadata
defval('metadatadir', '/Volumes/AOGquake/work/Ivy/CAATEX/metadata');
defval('datadir', '/Volumes/AOGquake/vault/CAATEX/')

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

    % determine length of time series data, time and space sampling
    [x, t, jdn] = HMrcvx(allfiles{ii});
    N = length(t);
    dt = t(2) - t(1);
    dz = hmdep(2) - hmdep(1);
    
    % construct 2D array storing signals from the hydrophone array
    xx = zeros([N length(snx)]);
    xx(:,1) = x';

    % loop over the hydrophones
    for jj = 2:length(snx)
        try
            x = HMrcvx(sprintf('%s/%03d/%s', masterdir, snx(jj), fname));
            xx(:,jj) = x';
        catch ME
            % if cannot read file, use data from the previous hydrophone
            ME.getReport()
            xx(:,jj) = xx(:,jj-1);
            continue
        end
    end

    % F-K transformation
    [XX, F, K] = fk(xx, dt, dz);

    % plot
    figure(1);
    set(gcf, 'Unit', 'inches', 'Position', [0 1 11 5]);
    gcf;
    ax1 = subplot('Position', [0.06 0.12 0.41 0.8]);
    imagesc(K, F, 10*log10(abs(XX).^2));
    grid on
    axis xy
    xlabel('wavenumber (m^{-1})')
    ylabel('frequency (Hz)')
    title(sprintf('station: %s - %s', station, datestr(jdn)))
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
    ylim([-50 50])
    cb = colorbar;
    cb.TickDirection = 'out';
    cb.Label.String = '10 log_{10}(spectral density) (dB)';
    set(gca, 'TickDir', 'out', 'Box', 'on', 'FontSize', 12)

    set(gcf, 'Renderer', 'painters')
    savename = sprintf('%s_%s_%s', mfilename, station, ...
        replace(fname, 'nc', 'eps'));
    figdisp(savename, [], [], 2, [], 'epstopdf')
end
end