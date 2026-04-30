function varargout = polarbasemap(ax, rlim, type)
% POLARBASEMAP(rlim)
% POLARBASEMAP(rlim, type)
% POLARBASEMAP(ax, rlim, ...)
% axhandl = POLARBASEMAP(...)
%
% Makes a basemap centered at the northpole. The basemap contains the land
% and the ocean as well as basic latitude/longitude grid. The coordinates
% of the plot remains cartesian with the northpole centered at (0,0) and
% the prime meridian (longitude=0) is pointing in -Y direction (downward).
% Use PLOTPOLARMAP to plot anything to the basemap. Otherwise, make sure 
% that you use [x,y] = (90-lat)*[sin(lon*pi/180), -cos(lon*pi/180]
% coordinates when you plot anything to the basemap.
%
% INPUT:
% ax        target axes                                     [default: gca]
% rlim      map radius (colatitude) limit in degrees        [default: 40]
% type      type of the map, either
%           'cartoon' -- land polygons on top of ocean
%           'topo'    -- GEBCO bathymetry map               [default]
%
% OUTPUT:
% axhandl   axes handle
%
% SEE ALSO:
% PLOTPOLARMAP
%
% Last modified by spipatprathanporn@ucsd.edu, 04/30/2026

defval('ax', gca)
defval('rlim', 40)
defval('type', 'topo')

% handle input parameters where the target axes is missing
if isnumeric(ax) && isscalar(ax)
    type = rlim;
    rlim = ax;
    ax = gca;
end

% list of longitudes for plotting the map rim and latitude grid
lon_rad = (0:0.025:2)' * pi;

hold(ax, 'on');

% cartoon mode - land polygons
if strcmp(type, 'cartoon')
    % background polygon = ocean
    pgon = polyshape(rlim*sin(lon_rad(1:end-1)), -rlim*cos(lon_rad(1:end-1)));
    plot(ax, pgon, 'FaceColor', [0.8 0.9 1], 'FaceAlpha', 1);
    
    % Maing coastlines

    % reading coastline file and plot on the tempolary axes. We don't need
    % the plot, just the coordinates
    fig_coast = figure;    
    [~, handl_coast] = plotcont([0 90], [360 90-rlim]);

    % convert (lon,lat) to (x,y)
    lon_coast_rad = deg2rad(handl_coast.XData');
    colat_coast = 90 - handl_coast.YData';
    delete(fig_coast);
    x_coast = colat_coast .* sin(lon_coast_rad);
    y_coast = -colat_coast .* cos(lon_coast_rad);
    xy_coast = [x_coast, y_coast];

    % split the coastlines at NaN entries
    lines_coast = splitnan(xy_coast);

    % Combine adjacent lines: they may overlap by one point or two points
    len = length(lines_coast);
    for ii = len:-1:2
        for jj = ii-1:-1:1
            % check if the sections overlap by 1 element
            if all(indeks(lines_coast{ii}, '1,:') == indeks(lines_coast{jj}, 'end,:'))
                lines_coast{jj} = [lines_coast{jj}; indeks(lines_coast{ii}, '2:end,:')];
                lines_coast(ii) = [];
                break
            elseif all(indeks(lines_coast{ii}, '1,:') == indeks(lines_coast{jj}, '1,:'))
                lines_coast{jj} = [flipud(lines_coast{jj}); indeks(lines_coast{ii}, '2:end,:')];
                lines_coast(ii) = [];
                break
            elseif all(indeks(lines_coast{ii}, 'end,:') == indeks(lines_coast{jj}, 'end,:'))
                lines_coast{jj} = [lines_coast{jj}; indeks(lines_coast{ii}, 'end-1:-1:1,:')];
                lines_coast(ii) = [];
                break
            elseif all(indeks(lines_coast{ii}, 'end,:') == indeks(lines_coast{jj}, '1,:'))
                lines_coast{jj} = [flipud(lines_coast{jj}); indeks(lines_coast{ii}, 'end-1:-1:1,:')];
                lines_coast(ii) = [];
                break
            end
            % check if the sections overlap by 2 elements
            if size(lines_coast{ii}, 1) < 2
                break
            end
            if all(indeks(lines_coast{ii}, '2,:') == indeks(lines_coast{jj}, 'end,:'))
                lines_coast{jj} = [lines_coast{jj}; indeks(lines_coast{ii}, '3:end,:')];
                lines_coast(ii) = [];
                break
            elseif all(indeks(lines_coast{ii}, '2,:') == indeks(lines_coast{jj}, '1,:'))
                lines_coast{jj} = [flipud(lines_coast{jj}); indeks(lines_coast{ii}, '3:end,:')];
                lines_coast(ii) = [];
                break
            elseif all(indeks(lines_coast{ii}, 'end-1,:') == indeks(lines_coast{jj}, 'end,:'))
                lines_coast{jj} = [lines_coast{jj}; indeks(lines_coast{ii}, 'end-2:-1:1,:')];
                lines_coast(ii) = [];
                break
            elseif all(indeks(lines_coast{ii}, 'end-1,:') == indeks(lines_coast{jj}, '1,:'))
                lines_coast{jj} = [flipud(lines_coast{jj}); indeks(lines_coast{ii}, 'end-2:-1:1,:')];
                lines_coast(ii) = [];
                break
            end
        end
    end
    
    % remove the point that is go pass a closed loop
    for ii = 1:length(lines_coast)
        if (size(lines_coast{ii}, 1) >= 2) && all(indeks(lines_coast{ii}, '2,:') == indeks(lines_coast{ii}, 'end,:'))
            lines_coast{ii} = indeks(lines_coast{ii}, '1:end-1,:');
        end
    end
    
    % fill-in the rim if the section ends at the outer rim
    for ii = 1:length(lines_coast)
        % check if the line is not a closed loop
        pt1 = indeks(lines_coast{ii}, '1,:');
        ptN = indeks(lines_coast{ii}, 'end,:');
        if any(pt1 ~= ptN)
            % check if both ends of the line are near the rim
            colat_pt1 = norm(pt1);
            colat_ptN = norm(ptN);
            if (rlim - colat_pt1 < 1) && (rlim - colat_ptN < 1)
                xrim = linspace(ptN(1), pt1(1), 61)';
                yrim = linspace(ptN(2), pt1(2), 61)';
    
                % rescale intermediate points to the rim
                drim = vecnorm([xrim, yrim], 2, 2);
                xrim = xrim * rlim ./ drim;
                yrim = yrim * rlim ./ drim;
    
                % fill-in the rim
                lines_coast{ii} = [lines_coast{ii}; [xrim, yrim]];
            end
        end
    end
    
    for ii = 1:length(lines_coast)
        coast = lines_coast{ii};
        if all(coast(1,:) == coast(end,:))
            pgon = polyshape(coast(1:end-1,1), coast(1:end-1,2));
            plot(ax, pgon, 'FaceColor', [0.7 0.9 0.6], 'FaceAlpha', 1)
        else
            pgon = polyshape(coast(:,1), coast(:,2));
            plot(ax, pgon, 'FaceColor', [0.7 0.9 0.6], 'FaceAlpha', 1)
        end
    end
    % Old code: uncomment this when you break anything after reading the
    % coastline file
    % plot(ax, colat_coast .* sin(lon_coast_rad), ...
    %     -colat_coast .* cos(lon_coast_rad), ...
    %     'Color', 'k', 'LineWidth', 1);

% plotting GEBCO bathymetry map
else
    x = linspace(-rlim, rlim, 1001);
    y = x;
    [xx, yy] = meshgrid(x, y);
    llat = 90 - sqrt(xx.^2 + yy.^2);
    llon = rad2deg(atan2(xx, -yy));

    % read the GEBCO file
    [lons, lats, elev] = bathymetry([], [-180 180], [90-rlim 90], false);
    topo = interp2(mod(lons+180,360)-180, lats, double(elev'), ...
        reshape(llon, [1 numel(llon)]), reshape(llat, [1 numel(llon)]), ...
        'linear');
    topo = reshape(topo, size(xx));
    
    im = imagesc([-rlim rlim], [-rlim rlim], topo);

    % set the colormap to suit the topography
    cb = cax2dem([-5000 5000]);
    delete(cb)

    % remove the color outside of the rim
    set(im, 'AlphaData', ~isnan(topo));
end

% plate boundaries
fig_plate = figure;
[~, XY] = plotplates([0 90], [360 90-rlim]);
delete(fig_plate);
lon_plate_rad = deg2rad(XY(:,1));
colat_plate = 90 - XY(:,2);
plot(ax, colat_plate .* sin(lon_plate_rad), ...
    -colat_plate .* cos(lon_plate_rad), ...
    'Color', 'r', 'LineWidth', 1, 'DisplayName', 'plate boundary');

% latitude grid
lat_grid = (10:10:rlim)';
if length(lat_grid) < 3
    lat_grid = (5:5:rlim)';
end
if length(lat_grid) < 3
    lat_grid = (1:1:rlim)';
end
for ii = 1:length(lat_grid)
    plot(ax, lat_grid(ii)*sin(lon_rad), -lat_grid(ii)*cos(lon_rad), ...
        'Color', [0.5 0.5 0.5], 'LineWidth', 0.5, 'HandleVisibility', 'off')
end

% longitude grid
lon_grid = (0:5) * pi/6;
for ii = 1:length(lon_grid)
    plot(ax, [-1 1] * rlim * sin(lon_grid(ii)), ...
        [1 -1] * rlim * cos(lon_grid(ii)), 'Color', [0.5 0.5 0.5], ...
        'LineWidth', 0.5, 'HandleVisibility', 'off')
end

% outer rim
plot(ax, rlim*sin(lon_rad), -rlim*cos(lon_rad), 'Color', [0 0 0], ...
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

% collecting output
varns = {ax};
varargout = varns(1:nargout);
end