function evs_str = readquerystack(ddir)
% evs_str = READQUERYSTACK(ddir)
%
% Reads all savefiles from QUERYSTACK in the directory and combine all of
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
% QUERYSTACK, STALTASTACK
%
% Last modified by spipatprathanporn@ucsd.edu, 10/12/2026

fnames = ls2cell(fullfile(ddir, 'querystack_output_*.mat'), 1)';

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

% remove duplicate events
dts = datetime(evs_str.PreferredTime);
[~, iu] = unique(dts);
for ii = 1:length(fn)
    evs_str.(fn{ii}) = evs_str.(fn{ii})(iu);
end
end