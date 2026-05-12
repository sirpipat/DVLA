function [evs_str, dd, outputdir] = querysot(lonlim, latlim, starttime, ...
    endtime, minmag, maxmag ,pstation, plocation, tstation, tlocation, name)
% [evs_str, dd, outputdir] = querysot([minlon maxlon], [minlon maxlon], ...
%     starttime, endtime, minmag, maxmag, pstation, [plon plat], ...
%     tstation, [tlon tlat], name)
%
% Query events inside [minlon maxlon] and [minlat maxlat], waveforms, and 
% response functions at pstation (P-wave) and tstation (T-wave).
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
%
% OUTPUT:
% evs_str               event information in a struct of arrays with a
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
%
% SEE ALSO:
% XCORRSOT
%
% Last modified by spipatprathanporn@ucsd.edu, 05/12/2025

outputdir = fullfile(getenv('IFILES'), 'SEISMOQUERY', name);

% base command for seismogram query
basecommand = sprintf('%s %s/irisFetch/seismo_query.py', ...
    getenv('PYTHON'), getenv('DVLA'));
basecommand = sprintf('%s --radius %f', basecommand, 5);

%% Part 1: query events
evs = irisFetch.Events('MinimumMagnitude', minmag, ...
    'MaximumMagnitude', maxmag, ...
    'MinimumLatitude', latlim(1), ...
    'MaximumLatitude', latlim(2), ...
    'MinimumLongitude', lonlim(1), ...
    'MaximumLongitude', lonlim(2), ...
    'StartTime', starttime, ...
    'EndTime', endtime);
N = length(evs);

% convert to struct of arrays
evs_str = array2struct(evs(end:-1:1));

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

%% Part 2: query the seismograms
for ii = 1:N
    dt_origin = datetime(evs_str.PreferredTime(ii), ...
        'Format', 'uuuu-MM-dd''T''HH:mm:ss.SSSSSS');
    dt_begin = dt_origin - minutes(1);
    dt_end = dt_origin + seconds(evs_str.tT(ii)) + minutes(10);

    % construct the Python call
    words = split(tstation, '.');
    if isempty(words{3})
        words{3} = '""';
    end
    command = sprintf('%s --lat %f --lon %f', basecommand, tlocation(2), tlocation(1));
    command = sprintf('%s --network %s --station %s', command, words{1}, words{2});
    command = sprintf('%s --start %s --end %s', command, dt_begin, dt_end);
    command = sprintf('%s --location %s', command, words{3}); 
    command = sprintf('%s --channels %s --outdir %s', command, ...
        words{4}, outputdir);
    command = sprintf('%s --format sac', command);

    % excecute the command
    if ~exist(fullfile(outputdir, sprintf('%s.%s_%s.sac', tstation, dt_begin, dt_end)), 'file') || ...
            ~exist(fullfile(outputdir, sprintf('%s.%s_%s.sacpz', tstation, dt_begin, dt_end)), 'file')
        system(command);
    else
        fprintf('%s.%s_%s.sac is already exist in %s\n', pstation, dt_begin, dt_end, outputdir);
    end

    % download local station seismograms for origin time correction
    words = split(pstation, '.');
    if isempty(words{3})
        words{3} = '""';
    end
    command = sprintf('%s --lat %f --lon %f', basecommand, plocation(2), plocation(1));
    command = sprintf('%s --network %s --station %s', command, words{1}, words{2});
    command = sprintf('%s --start %s --end %s', command, dt_begin, dt_end);
    command = sprintf('%s --location %s', command, words{3}); 
    command = sprintf('%s --channels %s --outdir %s', command, ...
        words{4}, outputdir);
    command = sprintf('%s --format sac', command);
    
    % excecute the command
    if ~exist(fullfile(outputdir, sprintf('%s.%s_%s.sac', pstation, dt_begin, dt_end)), 'file') || ...
            ~exist(fullfile(outputdir, sprintf('%s.%s_%s.sacpz', pstation, dt_begin, dt_end)), 'file')
        system(command);
    else
        fprintf('%s.%s_%s.sac is already exist in %s\n', pstation, dt_begin, dt_end, outputdir);
    end
end


end