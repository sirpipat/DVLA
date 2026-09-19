function varargout = querysot(lonlim, latlim, starttime, endtime, ...
    minmag, maxmag, pstation, plocation, tstation, tlocation, name)
% QUERYSOT([minlon maxlon], [minlat maxlat], starttime, endtime, ...
%     minmag, maxmag, pstation, [plon plat], tstation, [tlon tlat], ...
%     name)
% QUERYSOT(..., fname)
% QUERYSOT(..., evs_str_in)
% [evs_str, dd, outputdir, sname] = QUERYSOT(...)
%
% Query events inside [minlon maxlon] and [minlat maxlat], waveforms, and 
% response functions at pstation (P-wave) and tstation (T-wave).
% Alternatively, table of events from a file or as an event struct can be
% specified in variables "fname" and "evs_str_in", respectively. If
% neither of these is specified, it will call IRISFETCH.EVENTS to obtain
% the list of events. The outputs are also saved to a file. See the output
% variable sname for the name of the output file.
%
% INPUT:
% [minlon maxlon]       left and right boundary of the box
% [minlat maxlat]       top and bottom boundary of the box
% starttime             start time
% endtime               end time
% minmag                minimum magnitude
% maxmag                maximum magnitude
% pstation              network.station.location.channel for P-wave xcorr
% [plon plat]           lon,lat cooridnates of pstation
% tstation              network.station.location.channel for T-wave xcorr
% [tlon tlat]           lon,lat cooridnates of tstation
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
%       tP                          expected P-wave travel time to tstation
%       tT                          expected T-wave travel time to tstation
% dd                    event spatial separation in km
% outputdir             directory to the save files at: 
%                       $IFILES/SEISMOQUERY/name
% sname                 filename where the outputs and input arguments are
%                       saved: $IFILES/SEISMOQUERY/name/...
%                       querysot_output_[datetime("now")].mat
%
% SEE ALSO:
% READQUERYSOT, XCORRSOT
%
% Last modified by spipatprathanporn@ucsd.edu, 09/19/2026

% tracking the elapsed time
tic;

defval('fname', [])

outputdir = fullfile(getenv('IFILES'), 'SEISMOQUERY', name);
if ~exist(outputdir, 'dir')
    system(sprintf('mkdir %s', outputdir))
end

%% Part 1: query events
if ~isstruct(fname)
    if isempty(fname)
        evs = irisFetch.Events('MinimumMagnitude', minmag, ...
            'MaximumMagnitude', maxmag, ...
            'MinimumLatitude', latlim(1), ...
            'MaximumLatitude', latlim(2), ...
            'MinimumLongitude', lonlim(1), ...
            'MaximumLongitude', lonlim(2), ...
            'StartTime', starttime, ...
            'EndTime', endtime);
    else
        % read the event lists from the table
        T = readtable(fname);
        T = unique(T);
    
        % filter out out-of-range events
        starttime = datetime(starttime);
        endtime = datetime(endtime);
    
        wh = and(and(and(T.Latitude>=latlim(1), T.Latitude<=latlim(2)), ...
            and(T.Longitude>=lonlim(1), T.Longitude<=lonlim(2))), ...
            and(T.Date+T.Time>=starttime, T.Date+T.Time<=endtime));
        T = T(wh,:);
        
        % sort events
        T = sortrows(sortrows(T, 'Time', 'ascend'), 'Date', 'ascend');
        N = height(T);
    
        % place holder event list
        ev_nan = struct('PreferredTime', 'NaT', 'PreferredLatitude', NaN, ...
            'PreferredLongitude', NaN, 'PreferredDepth', NaN, ...
            'PreferredMagnitudeValue', NaN, 'PreferredMagnitudeType', '', ...
            'PublicId', 'evid=NaN');
        evs = repmat(ev_nan, [1 N]);
    
        for ii = 1:N
            ev = irisFetch.Events('MinmumMagnitude', T.Magnitude(ii)-0.5, ...
                'MaximumMagnitude', T.Magnitude(ii)+0.5, ...
                'MinimumLatitude', T.Latitude(ii)-0.5, ...
                'MaximumLatitude', T.Latitude(ii)+0.5, ...
                'MinimumLongitude', T.Longitude(ii)-0.5, ...
                'MaximumLongitude', T.Longitude(ii)+0.5, ...
                'StartTime', string(T.Date(ii)+T.Time(ii)-minutes(0.5), 'uuuu-MM-dd''T''HH:mm:ss.SSSSSS'), ...
                'EndTime', string(T.Date(ii)+T.Time(ii)+minutes(0.5), 'uuuu-MM-dd''T''HH:mm:ss.SSSSSS'));
        
            if length(ev)~=1
                keyboard
            else
                evs(ii) = ev;
            end
        end
    end
    
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
    evs_str = fname;
    N = length(evs_str.PreferredTime);
end

% compute the distance
dx = sin(deg2rad(mean(latlim))) .* ...
    deg2km(evs_str.PreferredLongitude - evs_str.PreferredLongitude');
dy = deg2km(evs_str.PreferredLatitude - evs_str.PreferredLatitude');
dz = evs_str.PreferredDepth - evs_str.PreferredDepth';
dd = sqrt(dx.^2 + dy.^2 + dz.^2);

% compute epicentral distance and the expected travel times
for ii = 1:N
    [evs_str.distkm(ii), evs_str.distdeg(ii)] = ...
        grcdist([evs_str.PreferredLongitude(ii) ...
        evs_str.PreferredLatitude(ii)], tlocation);
    tt = tauptime('mod', 'prem', 'dep', evs_str.PreferredDepth(ii), ...
        'ph', 'p,P', 'deg', evs_str.distdeg(ii));
    evs_str.tP(ii) = tt(1).time;
    evs_str.tT(ii) = evs_str.distkm(ii) / 1.51;
end

%% Part 2: constuct a query csv
query_str.network = repmat({''}, 2*N, 1);
query_str.station = repmat({''}, 2*N, 1);
query_str.location = repmat({''}, 2*N, 1);
query_str.channel = repmat({''}, 2*N, 1);
query_str.starttime = repmat({''}, 2*N, 1);
query_str.endtime = repmat({''}, 2*N, 1);
for ii = 1:N
    dt_origin = datetime(evs_str.PreferredTime(ii), ...
        'Format', 'uuuu-MM-dd''T''HH:mm:ss.SSSSSS');
    dt_begin = dt_origin - minutes(1);
    dt_end = dt_origin + seconds(evs_str.tT(ii)) + minutes(10);

    words = split(pstation, '.');
    % if isempty(words{3})
    %     words{3} = '""';
    % end
    query_str.network{2*ii-1} = words{1};
    query_str.station{2*ii-1} = words{2};
    query_str.location{2*ii-1} = words{3};
    query_str.channel{2*ii-1} = words{4};
    query_str.starttime{2*ii-1} = string(dt_begin);
    query_str.endtime{2*ii-1} = string(dt_end);

    words = split(tstation, '.');
    query_str.network{2*ii} = words{1};
    query_str.station{2*ii} = words{2};
    query_str.location{2*ii} = words{3};
    query_str.channel{2*ii} = words{4};
    query_str.starttime{2*ii} = string(dt_begin);
    query_str.endtime{2*ii} = string(dt_end);
end
query_T = struct2table(query_str);
csvname = fullfile(outputdir, 'querysot_query.csv');
writetable(query_T, csvname);

%% Part 3: query the seismograms
command = sprintf('%s %s/irisFetch/seismo_fetch.py', getenv('PYTHON'), ...
    getenv('DVLA'));
command = sprintf('%s --fname %s', command, csvname);
command = sprintf('%s --client IRIS --outdir %s', command, outputdir);
command = sprintf('%s --format sac', command);
system(command);

% for ii = 1:N
%     dt_origin = datetime(evs_str.PreferredTime(ii), ...
%         'Format', 'uuuu-MM-dd''T''HH:mm:ss.SSSSSS');
%     dt_begin = dt_origin - minutes(1);
%     dt_end = dt_origin + seconds(evs_str.tT(ii)) + minutes(10);
% 
%     % construct the Python call
%     words = split(tstation, '.');
%     if isempty(words{3})
%         words{3} = '""';
%     end
%     command = sprintf('%s --lat %f --lon %f', basecommand, tlocation(2), tlocation(1));
%     command = sprintf('%s --network %s --station %s', command, words{1}, words{2});
%     command = sprintf('%s --start %s --end %s', command, dt_begin, dt_end);
%     command = sprintf('%s --location %s', command, words{3}); 
%     command = sprintf('%s --channels %s --outdir %s', command, ...
%         words{4}, outputdir);
%     command = sprintf('%s --format sac', command);
% 
%     command = replace(command, '*', '\*');
%     command = replace(command, '?', '\?');
% 
%     % excecute the command
%     if ~exist(fullfile(outputdir, sprintf('%s.%s_%s.sac', tstation, dt_begin, dt_end)), 'file') || ...
%             ~exist(fullfile(outputdir, sprintf('%s.%s_%s.sacpz', tstation, dt_begin, dt_end)), 'file')
%         system(command);
%     else
%         fprintf('%s.%s_%s.sac is already exist in %s\n', tstation, dt_begin, dt_end, outputdir);
%     end
% 
%     % download local station seismograms for origin time correction
%     words = split(pstation, '.');
%     if isempty(words{3})
%         words{3} = '""';
%     end
%     command = sprintf('%s --lat %f --lon %f', basecommand, plocation(2), plocation(1));
%     command = sprintf('%s --network %s --station %s', command, words{1}, words{2});
%     command = sprintf('%s --start %s --end %s', command, dt_begin, dt_end);
%     command = sprintf('%s --location %s', command, words{3}); 
%     command = sprintf('%s --channels %s --outdir %s', command, ...
%         words{4}, outputdir);
%     command = sprintf('%s --format sac', command);
% 
%     command = replace(command, '*', '\*');
%     command = replace(command, '?', '\?');
% 
%     % excecute the command
%     if ~exist(fullfile(outputdir, sprintf('%s.%s_%s.sac', pstation, dt_begin, dt_end)), 'file') || ...
%             ~exist(fullfile(outputdir, sprintf('%s.%s_%s.sacpz', pstation, dt_begin, dt_end)), 'file')
%         system(command);
%     else
%         fprintf('%s.%s_%s.sac is already exist in %s\n', pstation, dt_begin, dt_end, outputdir);
%     end
% end

%% Part 3 save the output
sname = fullfile(outputdir, sprintf('%s_output_%s.mat', mfilename, ...
    datetime("now", "Format", "uuuu-MM-dd'T'HH:mm:ss.SSS")));
save(sname, 'evs_str', 'dd', 'outputdir', ...
    'lonlim', 'latlim', 'starttime', 'endtime', 'minmag', 'maxmag', ...
    'pstation', 'plocation', 'tstation', 'tlocation')
fprintf('The output is saved to %s\n', sname);

elapsed_time = toc;
fprintf('Elapsed time: %.2f seconds or %g events/minute\n', ...
    elapsed_time, N/(elapsed_time/60))

%% Part 4 collect the output
outputs = {evs_str, dd, outputdir, sname};
varargout = outputs(1:nargout);
end