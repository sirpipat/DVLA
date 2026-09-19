function evs_str = readquerysot(ddir)
% evs_str = READQUERYSOT(ddir)
%
% Reads all savefiles from QUERYSOT in the directory and combine all of
% them to a single struct.
%
% INPUT:
% ddir          directory to the QUERYSOT output files
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
%
% SEE ALSO:
% QUERYSOT, XCORRSOT
%
% Last modified by spipatprathanporn@ucsd.edu, 09/19/2026

fnames = ls2cell(fullfile(ddir, 'querysot_output_*.mat'), 1)';

load(fnames{1}, 'evs_str');
fn = fieldnames(evs_str);
for ii = 2:length(fnames)
    s = load(fnames{ii}, 'evs_str');
    for jj = 1:length(fn)
        if ~isfield(s.evs_str, fn{jj})
            s.evs_str.(fn{jj}) = repmat({''}, size(s.evs_str.(fn{1})));
        end
        evs_str.(fn{jj}) = horzcat(evs_str.(fn{jj}), s.evs_str.(fn{jj}));
    end
end

% make sure every datetime format is ended with ss.SSSSSS
N = length(evs_str.PreferredTime);
for ii = 1:N
    dt_str = evs_str.PreferredTime{ii};
    if length(dt_str) < 26
        if length(dt_str) == 19
            evs_str.PreferredTime{ii} = [dt_str '.000000'];
        else
            evs_str.PreferredTime{ii} = [dt_str  ...
                repmat('0', 1, 26-length(dt_str))];
        end
    end
end
end