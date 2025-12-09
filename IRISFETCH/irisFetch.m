classdef irisFetch

    
% IRISFETCH allows seamless access to data stored within the IRIS-DMC via FDSN services
%
% Spiritual successor of IRIS's IRISFETCH as the official IRISFETCH is
% no longer supported for Matlab 2023a onwards. In order to keep
% calling IRISFETCH in my codes, I need to make a replacement. Instead
% of using the IRIS Web Services Library java jar, it calls Python
% wrappers underhood. Make sure you install Python and Obspy package.
%
% This code is still under development and I do not guarantee that it
% will work in the same way as the IRIS's IRISFETCH.
%
% irisFetch Methods: only methods with asterisks (*) are implemented
%
% irisFetch waveform retrieval Methods:
%    Traces - retrieve sac-equivalent waveforms with channel metadata
%    SACfiles - as Traces above, but saves directly to a SAC file.
%
% irisFetch FDSN station webservice Methods:
%    Channels - retrieve metadata as an array of channels
%    Stations - retrieve metadata as an array of stations
%    Networks - retrieve metadata as an array of networks
%
% irisFetch FDSN event webservice Methods:
%  * Events - retrieve events parameters (such as origins and magnitudes) from a catalog
%
% irisFetch miscellaneous Methods:
%    Resp - retrive RESP formatted response data from the irisws-resp service
%    version - display the current version number
%    connectToJar - attempt to connect to the required IRIS-WS JAR file
%    runExamples - displays and runs some sample queries to the web service
%    Trace2SAC - writes a trace structure to a SAC file
%    SAC2Trace - reads one or more locally stored SAC files to a trace structure
%
%  irisFetch requires version 2.0 or greater of the IRIS Web Services Library java jar
%  for more details, click on 'connectToJar' above.
%
%  For additional guidance, type help <method>, use irisFetch.runExamples, or check out
%  the online manual https://ds.iris.edu/dms/nodes/dmc/software/downloads/irisFetch.m/
%
% Last modified by spipatprathanporn@ucsd.edu, 12/09/2025


methods(Static)
    function [events, urlParams] = Events(varargin)
        %irisFetch.Events retrieves event data from the IRIS-DMC
        %  ev = irisFetch.Events(param, value [, ...]) retrieves event data from the
        %  IRIS-DMC database as a matlab structure.  An arbitrary number of
        %  parameter-value pairs may be specified in order to narrow down the search
        %  results.
        %
        %  [ev, myParams] = irisFetch.Events( ... ) also returns the URL parameters that
        %  were used to make the query.
        %
        %  irisFetch.Events(..., 'BASEURL',alternateURL) specifies an alternate base URL
        %  to query. By default queries go to: http://service.iris.edu/fdsnws/event/1/
        %
        %  Usable parameters are listed below.  For detailed descriptions of their effect
        %  and use, consult the webservice webpage for events, available at:
        %
        %  http://service.iris.edu/fdsnws/event/1/
        %
        %  Examples:
        %
        %  Retrieve event parameters regarding earthquakes of a specific size and location:
        %      ev = irisFetch.Events('MinimumMagnitude',6.0,...
        %          'minimumLatitude',45,'maximumLatitude', 60,...
        %          'minimumLongitude', -150,'maximumLongitude', -90)
        %
        %PARAMETER LIST
        % contributor, endTime, eventId, fetchLimit, latitude, longitude, magnitudeType,
        % maximumDepth, maximumLatitude, maximumLongitude, maximumMagnitude,
        % minimumDepth, minimumLatitude, minimumLongitude, minimumMagnitude,
        % minimumRadius, maximumRadius, offset, startTime updatedAfter
        %
        %CONVENIENCE PARAMETERS
        %   'boxcoordinates'    : [minLat, maxLat, minLon, maxLon]   % use NaN as a wildcard
        %   'radialcoordinates' : [Lat, Lon, MaxRadius, MinRadius]   % MinRadius is optional
        %
        %NOTE: Any parameter-value pair that is not listed here will be passed along to
        %      the Event Service. this is to accomodate other datacenters that may have
        %      specialized request parameters.
        
        if nargin < 2
            warning('Need at least 1 parameter but there is none')
            events = [];
            urlParams = [];
            return
        end
        
        % construct the Python call
        command = sprintf('/Users/spipatprathanporn/anaconda3/bin/python %s/irisFetch/irisFetch_events.py', getenv('DVLA'));
        ii = 1;
        while ii < length(varargin)
            switch lower(varargin{ii})
                case 'starttime'
                    command = sprintf('%s --starttime %s', command, replace(varargin{ii+1}, ' ', 'T'));
                case 'endtime'
                    command = sprintf('%s --endtime %s', command,  replace(varargin{ii+1}, ' ', 'T'));
                case 'minimumlatitude'
                    command = sprintf('%s --minlatitude %s', command, string(varargin{ii+1}));
                case 'maximumlatitude'
                    command = sprintf('%s --maxlatitude %s', command, string(varargin{ii+1}));
                case 'minimumlongitude'
                    command = sprintf('%s --minlongitude %s', command, string(varargin{ii+1}));
                case 'maximumlongitude'
                    command = sprintf('%s --maxlongitude %s', command, string(varargin{ii+1}));
                case 'latitude'
                    command = sprintf('%s --latitude %s', command, string(varargin{ii+1}));
                case 'longitude'
                    command = sprintf('%s --longitude %s', command, string(varargin{ii+1}));
                case 'minimumradius'
                    command = sprintf('%s --minradius %s', command, string(varargin{ii+1}));
                case 'maximumradius'
                    command = sprintf('%s --maxradius %s', command, string(varargin{ii+1}));
                case 'minimumdepth'
                    command = sprintf('%s --mindepth %s', command, string(varargin{ii+1}));
                case 'maximumdepth'
                    command = sprintf('%s --maxdepth %s', command, string(varargin{ii+1}));
                case 'minimummagnitude'
                    command = sprintf('%s --minmagnitude %s', command, string(varargin{ii+1}));
                case 'maximummagnitude'
                    command = sprintf('%s --maxmagnitude %s', command, string(varargin{ii+1}));
                case 'eventid'
                    command = sprintf('%s --eventid %s', command, string(varargin{ii+1}));
                case 'boxcoordinates'
                    boxcoordinates = varargin{ii+1};
                    if ~isempty(boxcoordinates) && ~isnan(boxcoordinates(1))
                        command = sprintf('%s --minlatitude %s', command, string(boxcoordinates(1)));
                    end
                    if length(boxcoordinates) >= 2 && ~isnan(boxcoordinates(2))
                        command = sprintf('%s --maxlatitude %s', command, string(boxcoordinates(2)));
                    end
                    if length(boxcoordinates) >= 3 && ~isnan(boxcoordinates(3))
                        command = sprintf('%s --minlongitude %s', command, string(boxcoordinates(3)));
                    end
                    if length(boxcoordinates) >= 4 && ~isnan(boxcoordinates(4))
                        command = sprintf('%s --maxlongitude %s', command, string(boxcoordinates(4)));
                    end
                case 'radialcoordinates'
                    radialcoordinates = varargin{ii+1};
                    command = sprintf('%s --latitude %s', command, string(radialcoordinates(1)));
                    command = sprintf('%s --longitude %s', command, string(radialcoordinates(2)));
                    command = sprintf('%s --maxradius %s', command, string(radialcoordinates(3)));
                    if length(radialcoordinates) >= 4
                        command = sprintf('%s --minradius %s', command, string(radialcoordinates(4)));
                    end
            end
            ii = ii + 2;
        end
        outputfile = fullfile(getenv("IFILES"), 'MATFILES', 'event_catalog.mat');
        command = sprintf('%s --output %s', command, outputfile);

        % execute
        try
            [exitStatus, message] = system(command);
        catch ME
            ME.getReport;
        end

        % read the output file and structure
        events = [];
        urlParams = [];

        if exitStatus == 0
            load(outputfile, 'events');
            events = cell2mat(events);
            fprintf('%d Event(s) found\n', length(events));
            for ii = 1:length(events)
                fprintf('%s | %7.3f, %8.3f | %3.1f %s\n', ...
                    replace(events(ii).PreferredTime, ' ', 'T'), ...
                    events(ii).PreferredLatitude, ...
                    events(ii).PreferredLongitude, ...
                    events(ii).PreferredMagnitudeValue, ...
                    events(ii).PreferredMagnitudeType)
            end
        else
            fprintf(message);
            urlParams = message;
        end
    end
end

end