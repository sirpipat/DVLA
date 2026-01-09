function [t, x, jdn, num, depth] = getdvlaseis(station, hydrophone, dt_begin, dt_end, ...
    datadir, metadatadir)
% [t, x, jdn, num, depth] = GETDVLASEIS(station, hydrophone, dt_begin, dt_end)
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
% dt_begin      beginning datetime or datestr
% dt_begin      ending datetime or datestr
% 
% OUTPUT:
% t             time relative to jdn
% x             acoustic pressure
% jdn           datenum
% num           hydrophone number requested
% depth         depth of the hydrophone
% 
% Last modified by spipatprathanporn@ucsd.edu, 01/08/2026

% TODO: change to environment variables
% TODO: decide where to host data+metadata
defval('metadatadir', '/Volumes/aog/AOGquake/work/Ivy/CAATEX/metadata');
defval('datadir', '/Volumes/aog/AOGquake/vault/CAATEX/')

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
        datadir, metadatadir);
    return
elseif strcmpi(hydrophone, 'bottom')
    fprintf('Request seismogram from hydrophone number %d.\n', snx(end));
    [t, x, jdn, num, depth] = getdvlaseis(station, snx(end), dt_begin, dt_end, ...
        datadir, metadatadir);
    return
elseif strcmpi(hydrophone, 'middle')
    fprintf('Request seismogram from hydrophone number %d.\n', ...
        snx(ceil(length(snx)/2)));
    [t, x, jdn, num, depth] = getdvlaseis(station, snx(ceil(length(snx)/2)), ...
        dt_begin, dt_end, datadir, metadatadir);
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

% slice the seismogram
wh = and(dt >= dt_begin, dt <= dt_end);
x = x(wh);
t = t(wh);
num = hydrophone;
depth = hmdep(snx == hydrophone);
end