function [lons, lats, elev] = bathymetryglobal
% [lons, lats, elev] = BATHYMETRYGLOBAL
%
% Reads the downsampled GEBCO bathymetry grid for the global map, so it
% does not take forever. The resolution is 1/6 degrees.
%
% OUTPUT
% lons              longitudes of the requested bathymetry grid
% lats              latitudes of the requested bathymetry grid
% elev              elevations of the requested bathymetry grid
%                   Please note that elevation grid is Z = Z(lon,lat).
%
% SEE ALSO:
% BATHYMETRY
%
% Last modified by spipatprathanporn@ucsd.edu, 11/03/2026

fname = fullfile(getenv('IFILES'), 'TOPOGRAPHY', 'EARTH', ...
    'bathymetry_downsampled.mat');
load(fname, 'lons', 'lats', 'elev');
end