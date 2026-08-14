function varargout = polarbasemap2(ax, rlim, lon, lat, bath)
% POLARBASEMAP2(rlim)
% POLARBASEMAP2(rlim, lon, lat, bath)
% POLARBASEMAP2(ax, rlim, ...)
% axhandl = POLARBASEMAP2(...)
%
% Makes a basemap centered any given location. The basemap contains the 
% land and the ocean as well as basic distance/azimuth grid. The 
% coordinates of the plot remains cartesian centered at (0,0) and the 
% zeroth azimuth is pointing in Y direction (upward).
% Use PLOTPOLARMAP2 to plot anything to the basemap.
%
% INPUT:
% ax        target axes                                     [default: gca]
% rlim      map radius (colatitude) limit in degrees        [default: 40]
% lon       longitude of the center position                [default: 0]
% lat       latitude of the center position                 [default: 90]
% bath      whether to plot the bathymetry                  [default: true]
%
% OUTPUT:
% axhandl   axes handle
%
% SEE ALSO:
% PLOTPOLARMAP2, POLARBASEMAP
%
% Last modified by spipatprathanporn@ucsd.edu, 08/13/2026

defval('ax', gca)
defval('rlim', 40)
defval('lon', 0)
defval('lat', 90)
defval('bath', true)

% handle input parameters where the target axes is missing
if isnumeric(ax) && isscalar(ax)
    lat = lon;
    lon = rlim;
    rlim = ax;
    ax = gca;
end

cla

if lat == 90
    ax = polarbasemap(ax, rlim, 'topo');
else
    hold(ax, 'on');

    x = linspace(-rlim, rlim, 1001);
    y = x;
    [xx, yy] = meshgrid(x, y);
    rr = sqrt(xx.^2 + yy.^2);
    aa = mod(atan2(xx, yy) * 180/pi, 360);

    % flattens from 2D to 1D
    r = reshape(rr, numel(rr), []);
    az = reshape(aa, numel(aa), []);

    % only considers points within the map radius
    wh_in = (r <= rlim);
    r = r(wh_in);
    az = az(wh_in);

    % computes (lat,lon) of grid points on the map
    [lats2, lons2] = reckon(lat, lon, r, az);

    if abs(lat) + rlim < 90
        % rotates the center to be 0 meridian and sets the cutoff to be 180
        % so that the cutoff meridian will never be inside the map
        lons2_rotated = mod(lons2-lon+180, 360) - 180;
    
        % GEBCO query boundary
        lonleft = mod(min(lons2_rotated)+lon+180, 360) - 180;
        lonright = mod(max(lons2_rotated)+lon+180, 360) - 180;
    else
        lonleft = -180;
        lonright = 180;
    end
    latbottom = min(lats2);
    lattop = max(lats2);
    if bath
        % reads the bathymetry
        [lons, lats, elev] = bathymetry([], [lonleft lonright], [latbottom lattop]);
    
        % interpolate the bathymetry and make it 2D
        if lonleft < lonright
            topo = interp2(mod(lons+180,360)-180, lats, double(elev)', lons2, lats2);
        else
            topo = interp2(lons, lats, double(elev)', mod(lons2,360), lats2);
        end
        z = nan(size(wh_in));
        z(wh_in) = topo;
        z = reshape(z, sqrt(numel(rr)), []);
    
        % plots the elevation map
        im = imagesc([-rlim rlim], [-rlim rlim], z);
    
        % set the colormap to suit the topography
        cb = cax2dem([-8000 8000]);
        delete(cb)
    
        % remove the color outside of the rim
        set(im, 'AlphaData', rr <= rlim);
    else
        im = imagesc([-rlim rlim], [-rlim rlim], ones(size(rr)));
        colormap([1 1 1])

        % remove the color outside of the rim
        set(im, 'AlphaData', rr <= rlim);

        % coastlines
        fig_coast = figure;
        [~, ~, XYZ] = plotcont([min(mod(lons2,360)) lattop], [max(mod(lons2,360)) latbottom]);
        delete(fig_coast)
        wh = or(isnan(XYZ(:,1)), isnan(XYZ(:,2)));
        XYZ(wh,:) = nan;

        distdeg = distance(lat, lon, XYZ(:,2), XYZ(:,1), 'degrees');
        azrad = azimuth(lat, lon, XYZ(:,2), XYZ(:,1), 'degrees') * pi/180;
    
        distdeg(distdeg > rlim) = nan;
        azrad(distdeg > rlim) = nan;
    
        plot(ax, distdeg .* sin(azrad), distdeg .* cos(azrad), ...
            'Color', 'k', 'LineWidth', 1, 'DisplayName', 'coastlines');
    end

    % plate boundaries
    fig_plate = figure;
    [~, XY] = plotplates([min(mod(lons2,360)) lattop], [max(mod(lons2,360)) latbottom]);
    delete(fig_plate);
    % wh = find(and(XY(2:end,1)-XY(1:end-1,1)==0, XY(2:end,2)-XY(1:end-1,2)==0)) + 1;
    % for ii = length(wh):-1:1
    %     XY = [XY(1:wh(ii),:); nan nan; XY(wh(ii)+1:end,:)];
    % end
    wh = or(isnan(XY(:,1)), isnan(XY(:,2)));
    XY(wh,:) = nan;
    distdeg = distance(lat, lon, XY(:,2), XY(:,1), 'degrees');
    azrad = azimuth(lat, lon, XY(:,2), XY(:,1), 'degrees') * pi/180;

    distdeg(distdeg > rlim) = nan;
    azrad(distdeg > rlim) = nan;

    plot(ax, distdeg .* sin(azrad), distdeg .* cos(azrad), ...
        'Color', 'r', 'LineWidth', 1, 'DisplayName', 'plate boundary');

    % distance grid
    az_rad = (0:360)' * pi/180;
    dist_grid = (10:10:rlim)';
    if length(dist_grid) < 3
        dist_grid = (5:5:rlim)';
    end
    if length(dist_grid) < 3
        dist_grid = (1:1:rlim)';
    end
    for ii = 1:length(dist_grid)
        plot(ax, dist_grid(ii)*sin(az_rad), dist_grid(ii)*cos(az_rad), ...
            'Color', [0.5 0.5 0.5], 'LineWidth', 0.5, 'HandleVisibility', 'off')
    end
    
    % azimuth grid
    dist_grid = (0:5) * pi/6;
    for ii = 1:length(dist_grid)
        plot(ax, [-1 1] * rlim * sin(dist_grid(ii)), ...
            [1 -1] * rlim * cos(dist_grid(ii)), 'Color', [0.5 0.5 0.5], ...
            'LineWidth', 0.5, 'HandleVisibility', 'off')
    end
    
    % outer rim
    plot(ax, rlim*sin(az_rad), rlim*cos(az_rad), 'Color', [0 0 0], ...
        'LineWidth', 2, 'HandleVisibility', 'off');
    
    % axes data aspect ratio and sizing
    axis tight
    axis equal
    hold(ax, 'off')
    set(ax, 'FontSize', 12, 'Color', 'none', ...
        'XLim', [-1 1] * 1.02 * rlim, ...
        'YLim', [-1 1] * 1.02 * rlim)
    set(get(ax, 'XAxis'), 'Visible', 'off')
    set(get(ax, 'YAxis'), 'Visible', 'off')
end

% collecting output
set(ax, 'UserData', [lon lat])
varns = {ax};
varargout = varns(1:nargout);
end