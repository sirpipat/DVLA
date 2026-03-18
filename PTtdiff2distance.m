function gcarc = PTtdiff2distance(tdiff, evdp, bath, c, gcarc_T)
% gcarc = PTTDIFF2DISTANCE(tdiff, evdp, bath, c, gcarc_T)
%
% Estimates the epicentral distance from P-phase and T-phase travel time
% difference. It assumes the earthquakes is between 0 and 90 degrees away.
%
% INPUT:
% tdiff         (vector) of P-phase and T-phase travel time difference
% evdp          event (hypocenter) depth in km [Default: 10]
% bath          bathymetric depth at T-phase conversion point [Default: 0]
% c             oceaic acoustic wave speed in km/s [Default: 1.5]
% gcarc_T       distance on great-circle path in degrees from the epicenter
%               to P to T-phase conversion point [Default: 0]
%
% OUTPUT:
% gcarc         great-circle distance in degrees
%
% SEE ALSO:
% PTTRAVELTIME
%
% Last modified by spipatprathanporn@ucsd.edu, 03/17/2026

defval('evdp', 10)
defval('bath', 0)
defval('c', 1.5)
defval('gcarc_T', 0)

if size(tdiff, 2) > 1
    tdiff = tdiff';
end

x = (0:0.5:90)';
tPT = PTtraveltime(x, evdp, bath, c, gcarc_T);

gcarc = interp1(tPT, x, tdiff, 'linear');
end