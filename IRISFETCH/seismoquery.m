function varargout = seismoquery(event, radius, t_begin, t_end, ...
    channel, outputdir)
% [command, outputdir] = SEISMOQUERY(event, radius, t_begin, t_end, ...
%     channel, outputdir)
%
% Construct a command to query seimograms (SAC format) and their response 
% function (SACPZ format) from FDSN Web Service.
%
% INPUT:
% event         an event struct from IRISFETCH.EVENTS
%       - PreferredTime
%       - PreferredLatitude
%       - PreferredLongitude
%       - PreferredDepth
%       - PreferredMagnitudeValue
%       - PreferredMagnitudeType
%       - PublicId
% radius        maximum radius from the epicenter
% t_begin       beginning datetime  [default: 1 minute before origin time]
% t_end         ending datetime     [default: 1 hour after origin time]
% channel       channels to query   [default: 'BHZ']
% outputdir     where to save the SAC and SACPZ files
%
% OUTPUT:
% command       shell command for seismogram and response function query
% outputdir     where the SAC and SACPZ files are saved
%
% Last modified by spipatprathanporn@ucsd.edu, 03/25/2026

defval('radius', 20)
defval('channel', 'BHZ,BHN,BHE,BH1,BH2,BH3')
defval('outputdir', fullfile(getenv('IFILES'), 'SEISMOQUERY', ...
    sprintf('Event_%s', cindeks(split(event.PublicId, '='), 2))))

% calculate default starttime and endtime
origintime = datetime(event.PreferredTime, ...
    'Format', 'uuuu-MM-dd''T''HH:mm:ss.SSSSSS');
starttime = origintime - seconds(60);
endtime = origintime + hours(1);
defval('t_begin', starttime)
defval('t_end', endtime)

% make sure that the directory is absolute (i.e. no ~ for home directory)
outputdir = replace(outputdir, '~', getenv('HOME'));

% construct the Python call
command = sprintf('%s %s/irisFetch/seismo_query.py', getenv('PYTHON'), getenv('DVLA'));

command = sprintf('%s --lat %f --lon %f', command, ...
    event.PreferredLatitude, event.PreferredLongitude);
command = sprintf('%s --radius %f', command, radius);
command = sprintf('%s --start %s --end %s', command, t_begin, t_end);
command = sprintf('%s --channels %s --outdir %s', command, channel, ...
    outputdir);
command = sprintf('%s --format sac', command);

varns = {command, outputdir};
varargout = varns(1:nargout);
end