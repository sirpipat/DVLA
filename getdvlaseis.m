function [t, x, jdn, num, depth] = getdvlaseis(station, hydrophone, dt_begin, dt_end, fs, ...
    datadir, metadatadir)
% [t, x, jdn, num, depth] = GETDVLASEIS(station, hydrophone, dt_begin, dt_end, fs)
% 
% Obtains a section of seismogram from CAATEX hydrophone 
% array data set at a given station, hydrophone number, begin, and end time.
%
% INPUT:
% station       station name either 'NERSC1', 'NERSC2', 'NERSC3', 'SIO1', 
%               'SIO2', or 'SIO3'
% hydrophone    hydrophone equipment number, or these words
%     - 'list'      list valid hydrophone equipment number
%     - 'top'       get the top most hydrophone
%     - 'bottom'    get the bottom most hydrophone
%     - 'middle'    get the middle hydrophone
%     - 'all'       get all hydrophones
% dt_begin      beginning datetime or datestr
% dt_begin      ending datetime or datestr
% fs            sampling rate
% 
% OUTPUT:
% t             time relative to jdn
% x             acoustic pressure
% jdn           datenum
% num           hydrophone number requested
% depth         depth of the hydrophone
%
% EXAMPLES:
% % Example 1: list all valid hydrophone equipment number
% [nums, depths] = getdvlaseis('SIO1', 'list');
%
% % Example 2: get the time-series data from the top hydrophone
% [t, x, jdn, num, depth] = getdvlaseis('SIO1', nums(1), ...
%     '2019-12-15 12:00:00', datetime(2019, 12, 15, 12, 30, 0, 0));
% 
% Last modified by spipatprathanporn@ucsd.edu, 02/09/2026

% TODO: change to environment variables
% TODO: decide where to host data+metadata
defval('metadatadir', '/Volumes/aog/AOGquake/work/Ivy/CAATEX/metadata');
defval('datadir', '/Volumes/aog/AOGquake/vault/CAATEX/')
defval('fs', [])

% get metadata for the mooring
minfo = mooringinfo(station, metadatadir);
snx = minfo.snx;
hmdep = minfo.hmdep;

% determine which directories to read
masterdir = fullfile(datadir, station, 'HMrcv');

% verify equipment number
if strcmpi(hydrophone, 'list')
    fprintf('Return list of valid hydrophone numbers for station %s.\n', ...
        station);
    t = snx;
    x = hmdep;
    jdn = NaN;
    num = NaN;
    return
elseif strcmpi(hydrophone, 'top')
    fprintf('Request seismogram from hydrophone number %d.\n', snx(1));
    [t, x, jdn, num, depth] = getdvlaseis(station, snx(1), dt_begin, dt_end, ...
        fs, datadir, metadatadir);
    return
elseif strcmpi(hydrophone, 'bottom')
    fprintf('Request seismogram from hydrophone number %d.\n', snx(end));
    [t, x, jdn, num, depth] = getdvlaseis(station, snx(end), dt_begin, dt_end, ...
        fs, datadir, metadatadir);
    return
elseif strcmpi(hydrophone, 'middle')
    fprintf('Request seismogram from hydrophone number %d.\n', ...
        snx(ceil(length(snx)/2)));
    [t, x, jdn, num, depth] = getdvlaseis(station, snx(ceil(length(snx)/2)), ...
        dt_begin, dt_end, fs, datadir, metadatadir);
    return
elseif strcmpi(hydrophone, 'all')
    fprintf('Request seismograms from all hydrophones\n');
    if isempty(fs)
        fs = 50;
        fprintf('This option requires the sampling rate, fs\n');
        fprintf('It is set to the default value of 50 Hz\n');
    end
    % Requests the hydrophone from the first station to figure out the
    % length of the seismogram
    [t, x1, jdn] = getdvlaseis(station, snx(1), dt_begin, dt_end, ...
        fs, datadir, metadatadir);
    x = nan(length(snx), length(t));
    x(1,:) = x1;
    % Requests the rest of the mooring
    for ii = 2:length(snx)
        try
            [~, xi] = getdvlaseis(station, snx(ii), dt_begin, dt_end, ...
                fs, datadir, metadatadir);
            x(ii,:) = xi;
        catch ME
            ME.getReport()
            continue
        end
    end
    num = snx;
    depth = hmdep;
    return
end
if ~any(snx == hydrophone)
    error(['%d is not a valid hydrophone. The valid hydrophone number ' ...
        'for %s are %sand%d'], hydrophone, sprintf('%d ', snx(1:end-1)), ...
        snx(end));
end

% verify request time
dt_begin = datetime(dt_begin);
dt_end = datetime(dt_end);

% determine julian date
jdn = floor(days(dt_begin - datetime(2018,12,31,0,0,0)));
if dt_begin.Hour == 0 && dt_begin.Minute < 45
    jdn = jdn - 1;
    hours = 23;
else
    hours = 11;
end

% form the file name to read
fname = fullfile(masterdir, sprintf('%03d', hydrophone), ...
    sprintf('rcv_%d%d5950.nc', jdn, hours));
[x, t, jdn] = HMrcvx(fname);
dt = datetime(jdn, "ConvertFrom", "datenum", ...
    'Format', 'uuuu-MM-dd''T''HH:mm:ss.SSSSSS') + seconds(t);

if ~isempty(fs)
    % resample datetimes
    dt_sample = (dt_begin:seconds(1/fs):dt_end);

    fs0 = (length(t) - 1) / (t(end) - t(1));

    % lowpass filter to remove the potential alias
    x = lowpass(x, fs0, fs/2, 2, 1, 'butter', 'linear');

    % resample to the target sampling rate
    x = interp1(dt, x, dt_sample);
    
    % resample the times
    t = seconds(dt_sample - datetime(jdn, "ConvertFrom", "datenum", ...
        'Format', 'uuuu-MM-dd''T''HH:mm:ss.SSSSSS'));
    dt = dt_sample;
end

% slice the seismogram
wh = and(dt >= dt_begin, dt <= dt_end);
x = x(wh);
t = t(wh);
num = hydrophone;
depth = hmdep(snx == hydrophone);
end