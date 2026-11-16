function [lons, lats, isocean] = isoceanmap(reso, r, d, nolakes, ...
    nomediterranean)
% [lons, lats, isocean] = isoceanmap(reso, r, d, nolakes, nomediterranean)
%
% Returns the grid of points on the map with in distance (degrees) from the
% coast line with a given specified resolution (degrees).
%
% INPUT:
% reso                  grid resolution (degrees) [default: 1]
% r                     maximum distance (degrees) from the coastline 
%                       [default: 5]
% d                     ocean depth (m) [default: 0]
% nolakes               whether to remove lakes on the continents 
%                       [default: true]
% nomediterranean       whether to remove Mediterranean sea and Black sea
%                       [default: true]
%
% OUTPUT:
% lons              longitudes of the requested bathymetry grid
% lats              latitudes of the requested bathymetry grid
% elev              elevations of the requested bathymetry grid
%                   Please note that elevation grid is Z = Z(lat,lon).
%
% EXAMPLE:
% % simple call
% [lons, lats, isocean] = isoceanmap(1, 3, 100);
% % demo
% isoceanmap('demo');
%
% SEE ALSO:
% BATHYMETRYGLOBAL
%
% Last modified by spipatprathanporn@ucsd.edu, 11/16/2026

defval('reso', 1)
defval('r', 1)
defval('d', -0)
defval('nolakes', true)
defval('nomediterranean', true)

% demo
if strcmp(reso, 'demo')
    [lons, lats, elev] = bathymetryglobal;
    [lons2, lats2, isocean] = isoceanmap(1, 3, 100);

    figure(4)
    clf
    set(gcf, 'Units', 'inches', 'Position', [0 1 11 5.5])
    ax = subplot('Position', [0.08 0.10 0.88 0.85]);
    imagesc(lons, lats, elev');
    axis xy
    axis tight
    axis equal
    grid on
    cb = cax2dem([-1 1] * 8000);
    delete(cb)
    set(ax, 'FontSize', 11, 'Box', 'on', 'TickDir', 'out')
    xlabel('longitude')
    ylabel('latitude')
    title(sprintf(['%s demo | resolution: 1 degree | distance: ' ...
        '3 degrees | ocean depth: 100 m'], mfilename))
    xticks(0:30:360);
    yticks(-90:30:90);

    ax2 = addlayeraxes(ax);
    im = imagesc(lons2, lats2, isocean);
    axis xy
    axis tight
    axis equal
    grid on
    set(im, 'AlphaData', isocean * 0.5)
    colormap(ax2, [0 0 0; 0.5 0.5 1])
    set(ax2, 'FontSize', 11, 'Box', 'on', 'TickDir', 'out')
    set(gcf, 'Renderer', 'painters')
    exportgraphics(gcf, fullfile(getenv('EPS'), ...
        sprintf('%s_demo.pdf', mfilename)), 'ContentType', 'image', ...
        'Resolution', 600, 'Padding', 0);
    return
end

[lons, lats, elev] = bathymetryglobal;

if nolakes
    % removes lakes in Asia (Caspian Sea and Lake Baikal)
    elev(and(lons>43,lons<113),and(lats>31,lats<60)) = 10000;
    % removes lakes in North America
    % (The Great Lakes)
    elev(and(lons>266,lons<285),and(lats>42,lats<50)) = 10000;
    % (Basins in California)
    elev(and(lons>238,lons<239),and(lats>37.5,lats<38.5)) = 10000;
    elev(and(lons>243.5,lons<245),and(lats>32.5,lats<33.5)) = 10000;
    % removes lakes in South America
    elev(and(lons>287.5,lons<289),and(lats>8.5,lats<10.5)) = 10000;
    elev(and(lons>306,lons<309.75),and(lats>-2.5,lats<0.5)) = 10000;
    elev(and(lons>310,lons<312),and(lats>-2,lats<-1)) = 10000;
    % removes lakes in Australia
    elev(and(lons>136,lons<139),and(lats>-30,lats<-27)) = 10000;
    % removes lakes in North Africa
    elev(and(lons>20.5,lons<31),and(lats>28.5,lats<30.5)) = 10000;
    elev(and(lons>19,lons<20.5),and(lats>28.5,lats<30.1)) = 10000;
    elev(and(lons>5,lons<8.5),and(lats>33,lats<35)) = 10000;
end
if nomediterranean
    % removes anything in Mediterranean Sea and Black Sea
    elev(and(lons>0,lons<42),and(lats>30,lats<48)) = 10000;
end

% interpolates the grid to match the target resolution
lons2 = (reso:reso:360) - reso/2;
lats2 = (reso:reso:180) - reso/2 - 90;

[llons, llats]  = meshgrid(lons, lats);
[llons2,llats2] = meshgrid(lons2,lats2);

elev2 = interp2(llons, llats, double(elev'), llons2, llats2,'linear');

% constructs a circle window for minimum filter
n_pad = ceil(r / reso);
x = -n_pad:n_pad;
[yy, xx] = meshgrid(x, x);
w = ((xx*reso).^2 + (yy*reso).^2) > r^2;

% pads the elevation map before applying the minimum filter
elev3 = [elev2(:,end-n_pad+1:end) elev2 elev2(:,1:n_pad)];
elev3 = [repmat(1e5, [n_pad, size(elev3,2)]); elev3; ...
    repmat(1e5, [n_pad, size(elev3,2)])];

% applies the minimum filter
elev4 = nan(size(elev2));
for ii = 1:size(elev4,1)
    for jj = 1:size(elev4,2)
        minimap = elev3(ii:ii+2*n_pad, jj:jj+2*n_pad);
        minimap = minimap + 100000 * w;
        elev4(ii,jj) = min(minimap, [], 'all');
    end
end

% constructs the outputs
isocean = elev4 < -d;
lons = lons2;
lats = lats2;
end