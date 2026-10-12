function varargout = querystack(station, location, minradius, ...
    maxradius, minmag, maxmag, starttime, endtime, maxduration, name, ...
    evs_str_in)
% QUERYSTACK(tstation, [stlo stla], minradius, maxradius, minmag, ...
%     maxmag, starttime, endtime, maxduration, name)
% QUERYSTACK(..., fname)
% QUERYSTACK(..., evs_str_in)
% [evs_str, outputdir, sname] = QUERYSTACK(...)
%
% Query events inside [minradius maxradius] from tstation, waveforms, and 
% response functions at tstation (T-wave).
% Alternatively, table of events from a file or as an event struct can be
% specified in variables "fname" and "evs_str_in", respectively. If
% neither of these is specified, it will call IRISFETCH.EVENTS to obtain
% the list of events. The outputs are also saved to a file. See the output
% variable sname for the name of the output file.
%
% INPUT:
% tstation              network.station.location.channel for T-wave
% [stlo stla]           lon,lat cooridnates of tstation
% minradius             minimum radius
% maxradius             maximum radius
% minmag                minimum magnitude
% maxmag                maximum magnitude
% starttime             start time
% endtime               end time
% maxduration           maximum waveform duration in minutes
% name                  name you want to give for this query
% fname                 (optional) filename of a table of events
% evs_str_in            (optional) input event information in a struct of 
%                       arrays offollowing variables
%       PreferredTime               event time     as from irisFetch.Events
%       PreferredLatitude           event latitude
%       PreferredLongitude          event longitude      
%       PreferredDepth              event depth
%       PreferredMagnitudeValue     event magnitude value
%       PreferredMagnitudeType      event magnitude type
%
% OUTPUT:
% evs_str               event information in a struct of arrays with 
%                       following variables:
%       PreferredTime               event time     as from irisFetch.Events
%       PreferredLatitude           event latitude
%       PreferredLongitude          event longitude      
%       PreferredDepth              event depth
%       PreferredMagnitudeValue     event magnitude value
%       PreferredMagnitudeType      event magnitude type
%       PublicID                    event public ID
%       distkm                      distance to tstation in km
%       distdeg                     distance to tstation in degrees
%       tT                          expected T-wave travel time to tstation
% outputdir             directory to the save files at: 
%                       $IFILES/SEISMOQUERY/name
% sname                 filename where the outputs and input arguments are
%                       saved: $IFILES/SEISMOQUERY/name/...
%                       querysot_output_[datetime("now")].mat
%
% SEE ALSO:
% QUERYSTACK_ROUTINE, READQUERYSTACK, QUERYSOT
%
% Last modified by spipatprathanporn@ucsd.edu, 10/12/2026

% tracking the elapsed time
tic;

defval('fname', [])
defval('evs_str_in', [])

outputdir = fullfile(getenv('IFILES'), 'SEISMOQUERY', name);
if ~exist(outputdir, 'dir')
    system(sprintf('mkdir %s', outputdir))
end

%% Part 1: Query events
if isempty(evs_str_in)
    evs = irisFetch.Events('RadialCoordinates', ...
        [location(2) location(1) maxradius minradius], ...
        'MinimumMagnitude', minmag, ...
        'MaximumMagnitude', maxmag, ...
        'StartTime', starttime, ...
        'EndTime', endtime);
    
    N = length(evs);
        
    % sort events by time
    dt_origins = repmat(datetime('now','Format',...
        'uuuu-MM-dd''T''HH:mm:ss.SSSSSS'), N, 1);
    for ii = 1:N
        dt_origins(ii) = datetime(evs(ii).PreferredTime, ...
            'Format', 'uuuu-MM-dd''T''HH:mm:ss.SSSSSS');
    end
    [~,ii_sort] = sort(dt_origins);
    evs = evs(ii_sort);
    
    % convert to struct of arrays
    evs_str = array2struct(evs);
else
    evs_str = evs_str_in;
    N = length(evs_str.PreferredTime);
end

% compute epicentral distance and the expected travel times
for ii = 1:N
    [evs_str.distkm(ii), evs_str.distdeg(ii)] = ...
        grcdist([evs_str.PreferredLongitude(ii) ...
        evs_str.PreferredLatitude(ii)], location);
    evs_str.tT(ii) = evs_str.distkm(ii) / 1.51;
end

%% Part 2: constuct a query csv
query_str.network = repmat({''}, N, 1);
query_str.station = repmat({''}, N, 1);
query_str.location = repmat({''}, N, 1);
query_str.channel = repmat({''}, N, 1);
query_str.starttime = repmat({''}, N, 1);
query_str.endtime = repmat({''}, N, 1);
for ii = 1:N
    dt_origin = datetime(evs_str.PreferredTime(ii), ...
        'Format', 'uuuu-MM-dd''T''HH:mm:ss.SSSSSS');
    dt_begin = dt_origin - minutes(1);
    dt_end = dt_origin + minutes(maxduration);

    words = split(station, '.');
    % if isempty(words{3})
    %     words{3} = '""';
    % end
    query_str.network{ii} = words{1};
    query_str.station{ii} = words{2};
    query_str.location{ii} = words{3};
    query_str.channel{ii} = words{4};
    query_str.starttime{ii} = string(dt_begin);
    query_str.endtime{ii} = string(dt_end);
end
query_T = struct2table(query_str);
csvname = fullfile(outputdir, sprintf('%s_query.csv', mfilename));
writetable(query_T, csvname);

%% Part 3: query the seismograms
command = sprintf('%s %s/irisFetch/seismo_fetch.py', getenv('PYTHON'), ...
    getenv('DVLA'));
command = sprintf('%s --fname %s', command, csvname);
command = sprintf('%s --client IRIS --outdir %s', command, outputdir);
command = sprintf('%s --format sac', command);
system(command);

%% Part 4 save the output
sname = fullfile(outputdir, sprintf('%s_output_%s.mat', mfilename, ...
    datetime("now", "Format", "uuuu-MM-dd'T'HH:mm:ss.SSS")));
save(sname, 'evs_str', 'outputdir', 'station', 'location', ...
    'minradius', 'maxradius', 'minmag', 'maxmag', 'starttime', ...
    'endtime', 'maxduration', 'name', 'evs_str_in')
fprintf('The output is saved to %s\n', sname);

elapsed_time = toc;
fprintf('Elapsed time: %.2f seconds or %g events/minute\n', ...
    elapsed_time, N/(elapsed_time/60))

%% Part 5 collect the output
outputs = {evs_str, outputdir, sname};
varargout = outputs(1:nargout);
end