function [tdiff, tP, tT] = PTtraveltime(gcarc, evdp, bath, c, gcarc_T)
% [tdiff, tP, tT] = PTTRAVELTIME(gcarc, evdp, bath, c, gcarc_T)
% 
% Computes approximated P-phase and T-phase travel time difference. It
% assumes the station on the coast at zero elevation.
%
% INPUT:
% gcarc         (vector) of great-circle distance in degrees
% evdp          event (hypocenter) depth in km [Default: 10]
% bath          bathymetric depth at T-phase conversion point [Default: 0]
% c             oceaic acoustic wave speed in km/s [Default: 1.5]
% gcarc_T       distance on great-circle path in degrees from the epicenter
%               to P to T-phase conversion point [Default: 0]
%
% OUTPUT:
% tdiff         P-phase and T-phase travel time difference
% tP            P-phase travel time
% tT            T-phase travel time
%
% Last modified by spipatprathanporn@ucsd.edu, 03/17/2026

defval('evdp', 10)
defval('bath', 0)
defval('c', 1.5)
defval('gcarc_T', 0)


tP = nan(size(gcarc));
tT = nan(size(gcarc));

for ii = 1:size(gcarc,1)
        try
            tP(ii) = getfield(indeks(tauptime('dep',evdp,...
                'ph','p,P','deg',gcarc(ii),'mod','ak135'),1),'time');
        catch
        end
        % tT = p(hypocenter->conversion pt) + T(conversion pt->station)
        try
            tp = getfield(indeks(tauptime('dep',evdp,'stdp',bath,...
                'ph','p,P','deg',gcarc_T,'mod','ak135'),1),'time');
            tT(ii) = tp + deg2km(abs(gcarc(ii)-gcarc_T)) / c;
        catch
        end
end

tdiff = tT - tP;
end