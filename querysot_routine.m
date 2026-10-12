function varargout = querysot_routine(lonlim, latlim, starttime, ...
    endtime, minmag, maxmag, pstation, plocation, tstation, tlocation, ...
    name, fname, firstcall)
% QUERYSOT_ROUTINE([minlon maxlon], [minlat maxlat], starttime, ...
%     endtime, minmag, maxmag, pstation, [plon plat], tstation, ...
%     [tlon tlat], name)
% QUERYSOT_ROUTINE(..., fname)
% QUERYSOT_ROUTINE(..., evs_str_in)
% [evs_str, outputdir] = QUERYSOT_ROUTINE(...)
%
% Query events inside [minlon maxlon] and [minlat maxlat], waveforms, and 
% response functions at pstation (P-wave) and tstation (T-wave).
% Alternatively, table of events from a file or as an event struct can be
% specified in variables "fname" and "evs_str_in", respectively. If
% neither of these is specified, it will call IRISFETCH.EVENTS to obtain
% the list of events. In case the list of events is too large to query at 
% once, it will automatically divided into smaller chunks and query one at
% a time. The outputs are also saved to a file.
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
% outputdir             directory to the save files at: 
%                       $IFILES/SEISMOQUERY/name
%
% SEE ALSO:
% QUERYSOT, READQUERYSOT, XCORRSOT
%
% Last modified by spipatprathanporn@ucsd.edu, 10/12/2026

% tracking the elapsed time
timerval = tic;

defval('fname', [])
defval('firstcall', true)

outputdir = fullfile(getenv('IFILES'), 'SEISMOQUERY', name);
if ~exist(outputdir, 'dir')
    system(sprintf('mkdir %s', outputdir));
end

dt_start = datetime(starttime, 'Format', 'uuuu-MM-dd''T''HH:mm:ss');
dt_end = datetime(endtime, 'Format', 'uuuu-MM-dd''T''HH:mm:ss');

%% Schedule the query
try
    querysot(lonlim, latlim, starttime, endtime, minmag, maxmag, ...
        pstation, plocation, tstation, tlocation, name, fname);
catch ME
    if contains(ME.message, 'Index exceeds array bounds.')
        fprintf('Split task into half\n')
        dt_half = dt_start + 0.5 * (dt_end - dt_start);
        % don't immediately spam requests
        pause(10)
        % request the first half
        fprintf('Query events from %s to %s ...\n', starttime, ...
            string(dt_half));
        querysot_routine(lonlim, latlim, starttime, string(dt_half), ...
            minmag, maxmag, pstation, plocation, tstation, tlocation, ...
            name, fname, false);
        % don't immediately spam requests
        pause(10)
        % request the second half
        fprintf('Query events from %s to %s ...\n', string(dt_half), ...
            endtime);
        querysot_routine(lonlim, latlim, string(dt_half), endtime, ...
            minmag, maxmag, pstation, plocation, tstation, tlocation, ...
            name, fname, false);
    end
end

%% Part 5 collect the output
if nargout < 1
    return
end
evs_str = readquerysot(outputdir);
outputs = {evs_str, outputdir};
varargout = outputs(1:nargout);

if firstcall
    elapsed_time = toc(timerval);
    fprintf('%s Elapsed time: %.2f seconds or %g events/minute\n', ...
        mfilename, elapsed_time, ...
        length(evs_str.PreferredTime)/(elapsed_time/60))
end
end