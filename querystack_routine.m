function varargout = querystack_routine(station, location, minradius, ...
    maxradius, minmag, maxmag, starttime, endtime, maxduration, name, ...
    evs_str_in)
% QUERYSTACK_ROUTINE(tstation, [stlo stla], minradius, maxradius, ...
%     minmag, maxmag, starttime, endtime, maxduration, name)
% QUERYSTACK_ROUTINE(..., fname)
% QUERYSTACK_ROUTINE(..., evs_str_in)
% [evs_str, outputdir] = QUERYSTACK_ROUTINE(...)
%
% Routinely query events inside [minradius maxradius] from tstation, 
% waveforms, and response functions at tstation (T-wave). When the query is
% too much, or it times out, it divided the request to smaller chunks.
% Then, query one chunk at a time.
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
%
% SEE ALSO:
% QUERYSTACK, READQUERYSTACK, QUERYSOT_ROUTINE
%
% Last modified by spipatprathanporn@ucsd.edu, 10/12/2026

defval('evs_str_in', [])

outputdir = fullfile(getenv('IFILES'), 'SEISMOQUERY', name);
if ~exist(outputdir, 'dir')
    system(sprintf('mkdir %s', outputdir));
end

dt_start = datetime(starttime, 'Format', 'uuuu-MM-dd''T''HH:mm:ss');
dt_end = datetime(endtime, 'Format', 'uuuu-MM-dd''T''HH:mm:ss');

%% Schedule the query
try
    [~, outputdir] = querystack(station, location, ...
        minradius, maxradius, minmag, maxmag, starttime, endtime, ...
        maxduration, name, evs_str_in);
catch ME
    if contains(ME.message, 'Index exceeds array bounds.')
        fprintf('Split task into half\n')
        dt_half = dt_start + 0.5 * (dt_end - dt_start);
        % don't immediately spam requests
        pause(10)
        % request the first half
        fprintf('Query events from %s to %s ...\n', starttime, ...
            string(dt_half));
        querystack_routine(station, location, minradius, maxradius, ...
            minmag, maxmag, starttime, string(dt_half), maxduration, ...
            name, evs_str_in);
        % don't immediately spam requests
        pause(10)
        % request the second half
        fprintf('Query events from %s to %s ...\n', string(dt_half), ...
            endtime);
        querystack_routine(station, location, minradius, maxradius, ...
            minmag, maxmag, string(dt_half), endtime, maxduration, ...
            name, evs_str_in);
    end
end

%% Part 5 collect the output
if nargout < 1
    return
end
evs_str = readquerystack(outputdir);
outputs = {evs_str, outputdir};
varargout = outputs(1:nargout);
end