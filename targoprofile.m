function [x, z, T, dT] = targoprofile(lonlat1, lonlat2, npts, dt, fname)
% [x, z, T, dT] = TARGOPROFILE([lon1 lat1], [lon2 lat2], npts, dt, fname)
%
% Makes a cross section of ARGO mean temperature and temperature anomaly 
% along the path at a given datetime.
%
% INPUT:
% [lon1 lat1]       longitude and latitude of the starting point (degrees)
% [lon2 lat2]       longitude and latitude of the ending point (degrees)
% npts              number of points on the great-circle path
% dt                datetime
% fname             filename of the ARGO temperature file
%                   
% OUTPUT:
% x                 distance along the great circle path (km)
% z                 water depth (m)
% T                 mean temperature:     T = T(z, x)
% dT                temperature anomaly: dT = dT(z, x)
%
% Last modified by spipatprathanporn@ucsd.edu, 07/31/2026

% Please note that the longitude grid is from 20.5 to 379.5 degrees (cutoff
% at 20 degrees)
LON_CUTOFF = 20;

% default file name
defval('fname', fullfile(getenv('IFILES'), 'EARTHMODELS', 'PHYSICAL', ...
    'ARGO', 'RG_ArgoClim_Temperature_2019.nc'))

% compute the trajectory
[azDeg, ~, distDeg] = azimdist(lonlat1, lonlat2);
arclens = linspace(0, distDeg, npts);
[lats, lons] = reckon(lonlat1(2), lonlat1(1), arclens, azDeg);
lons = mod(lons - LON_CUTOFF, 360) + LON_CUTOFF; % 20 < lons < 380
lonlat1(1) = mod(lonlat1(1) - LON_CUTOFF, 360) + LON_CUTOFF;
lonlat2(1) = mod(lonlat2(1) - LON_CUTOFF, 360) + LON_CUTOFF;


% This handles cases where lons span across 20 E longitude
% Determine if there is any jump in lons
difflons = lons(2:end) - lons(1:(end-1));
if and(~or(all(difflons < 0), all(difflons > 0)), ~all(difflons == 0))
    is_across20 = true;
else
    is_across20 = false;
end


% Reads the list of longitudes, latitudes, pressures, and times of the
% whole temperature grids
LONGITUDES = double(ncread(fname, 'LONGITUDE')');
LATITUDES  = double(ncread(fname, 'LATITUDE')');
PRESSURES = double(ncread(fname, 'PRESSURE')');
TIMES = double(ncread(fname, 'TIME')');
DTIMES = datetime('16-JAN-2004 00:00:00') + calmonths(TIMES - 0.5);

% determines which elements to read
if is_across20
    lon1 = ceil(min(lons(1), lons(end)) + 0.5) - 0.5;
    lon2 = floor(max(lons(1), lons(end)) + 0.5) - 0.5;
    lat1 = fllor(min(lonlat1(2), lonlat2(2)) + 0.5) - 0.5;
    lat2 = ceil(max(lonlat1(2), lonlat2(2)) + 0.5) - 0.5;

    [~, ilon1] = min(abs(LONGITUDES - lon1));
    [~, ilon2] = min(abs(LONGITUDES - lon2));
    [~, ilat1] = min(abs(LATITUDES - lat1));
    [~, ilat2] = min(abs(LATITUDES - lat2));

    % reads the mean temperature (left half and right half)
    ARGO_TEMPERATURE_MEAN_LEFT = ncread(fname, 'ARGO_TEMPERATURE_MEAN', ...
        [ilon2 ilat1 1], [Inf ilat2-ilat1+1 Inf]);
    ARGO_TEMPERATURE_MEAN_RIGHT = ncread(fname, 'ARGO_TEMPERATURE_MEAN', ...
        [1 ilat1 1], [ilon1 ilat2-ilat1+1 Inf]);
    ARGO_TEMPERATURE_MEAN = cat(1, ARGO_TEMPERATURE_MEAN_LEFT, ...
        ARGO_TEMPERATURE_MEAN_RIGHT);
    
    % reads the temperature anomaly (left half and right half)
    ARGO_TEMPERATURE_ANOMALY_LEFT = ncread(fname, 'ARGO_TEMPERATURE_ANOMALY', ...
        [ilon2 ilat1 1 1], [Inf ilat2-ilat1+1 Inf Inf]);
    ARGO_TEMPERATURE_ANOMALY_RIGHT = ncread(fname, 'ARGO_TEMPERATURE_ANOMALY', ...
        [1 ilat1 1 1], [ilon1 ilat2-ilat1+1 Inf Inf]);
    ARGO_TEMPERATURE_ANOMALY = cat(1, ARGO_TEMPERATURE_ANOMALY_LEFT, ...
        ARGO_TEMPERATURE_ANOMALY_RIGHT);

    % longitudes and latitudes in the read temperature grids
    LON_USED = [LONGITUDES(ilon2:end) LONTIDUTES(1:ilon1)];
    LAT_USED = LATITUDES(ilat1:ilat2);
else
    lon1 = floor(min(lonlat1(1), lonlat2(1)) + 0.5) - 0.5;
    lon2 = ceil(max(lonlat1(1), lonlat2(1)) + 0.5) - 0.5;
    lat1 = floor(min(lonlat1(2), lonlat2(2)) + 0.5) - 0.5;
    lat2 = ceil(max(lonlat1(2), lonlat2(2)) + 0.5) - 0.5;

    [~, ilon1] = min(abs(LONGITUDES - lon1));
    [~, ilon2] = min(abs(LONGITUDES - lon2));
    [~, ilat1] = min(abs(LATITUDES - lat1));
    [~, ilat2] = min(abs(LATITUDES - lat2));

    % read the mean temperature
    ARGO_TEMPERATURE_MEAN = ncread(fname, 'ARGO_TEMPERATURE_MEAN', ...
        [ilon1 ilat1 1], [ilon2-ilon1+1 ilat2-ilat1+1 Inf]);

    % read the temperature anomaly
    ARGO_TEMPERATURE_ANOMALY = ncread(fname, 'ARGO_TEMPERATURE_ANOMALY', ...
        [ilon1 ilat1 1 1], [ilon2-ilon1+1 ilat2-ilat1+1 Inf Inf]);

    % longitudes and latitudes in the read temperature grids
    LON_USED = LONGITUDES(ilon1:ilon2);
    LAT_USED = LATITUDES(ilat1:ilat2);
end

% slice everything in time first
DTIME1 = max(DTIMES(DTIMES <= dt));
DTIME2 = min(DTIMES(DTIMES >= dt));
[~, it1] = min(abs(DTIMES - DTIME1));
[~, it2] = min(abs(DTIMES - DTIME2));
% linear interpolation
W1 = (dt - DTIME2) / (DTIME1 - DTIME2);
W2 = (dt - DTIME1) / (DTIME2 - DTIME1);
ARGO_TEMPERATURE_ANOMALY = W1 * ARGO_TEMPERATURE_ANOMALY(:,:,:,it1) + ...
    W2 * ARGO_TEMPERATURE_ANOMALY(:,:,:,it2);
ARGO_TEMPERATURE_ANOMALY = squeeze(ARGO_TEMPERATURE_ANOMALY);

% slice along the path
T = interpn(LON_USED, LAT_USED, PRESSURES, ARGO_TEMPERATURE_MEAN, ...
    repmat(lons, [1 length(PRESSURES)]), ...
    repmat(lats, [1 length(PRESSURES)]), repelem(PRESSURES, npts));
T = reshape(T, [npts length(PRESSURES)])';

dT = interpn(LON_USED, LAT_USED, PRESSURES, ARGO_TEMPERATURE_ANOMALY, ...
    repmat(lons, [1 length(PRESSURES)]), ...
    repmat(lats, [1 length(PRESSURES)]), repelem(PRESSURES, npts));
dT = reshape(dT, [npts length(PRESSURES)])';

% computes the vectors of x and z
x = deg2km(arclens);
% Converts pressures from dbar to deph in meters
% Assumes g=9.8 m/s^2, and saltwater density of 1020 kg/m^3
z = PRESSURES * 10000 / 9.8 / 1020;
end