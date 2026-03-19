function varargout = polarbasemap(ax, rlim)
% POLARBASEMAP(rlim)
% POLARBASEMAP(ax, rlim)
% axhandl = POLARBASEMAP(...)
%
% Makes a basemap centered at the northpole.
%
% INPUT:
% ax        target axes [default: gca]
% rlim      radius (colatitude) limit in degrees [default: 40]
%
% OUTPUT:
% axhandl   axes handle
%
% Last modified by spipatprathanporn@ucsd.edu, 03/19/2026

defval('ax', gca)
defval('rlim', 40)

if isnumeric(ax) && isscalar(ax)
    rlim = ax;
    ax = gca;
end

% background polygon
lon_rad = (0:0.025:2)' * pi;
pgon = polyshape(rlim*sin(lon_rad(1:end-1)), -rlim*cos(lon_rad(1:end-1)));
hold(ax, 'on');
plot(ax, pgon, 'FaceColor', [1 1 1], 'FaceAlpha', 1);

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
        'Color', [0.75 0.75 0.75], 'LineWidth', 0.5)
end

% longitude grid
lon_grid = (0:5) * pi/6;
for ii = 1:length(lon_grid)
    plot(ax, [-1 1] * rlim * sin(lon_grid(ii)), ...
        [1 -1] * rlim * cos(lon_grid(ii)), 'Color', [0.75 0.75 0.75], ...
        'LineWidth', 0.5)
end

% coastline
fig_coast = figure;
[~, handl_coast] = plotcont([0 90], [360 90-rlim]);
lon_coast_rad = deg2rad(handl_coast.XData);
colat_coast = 90 - handl_coast.YData;
delete(fig_coast);
plot(ax, colat_coast .* sin(lon_coast_rad), ...
    -colat_coast .* cos(lon_coast_rad), ...
    'Color', 'k', 'LineWidth', 1);

% plate boundaries
fig_plate = figure;
[~, XY] = plotplates([0 90], [360 90-rlim]);
delete(fig_plate);
lon_plate_rad = deg2rad(XY(:,1));
colat_plate = 90 - XY(:,2);
plot(ax, colat_plate .* sin(lon_plate_rad), ...
    -colat_plate .* cos(lon_plate_rad), ...
    'Color', 'r', 'LineWidth', 1);

% outer rim
plot(ax, rlim*sin(lon_rad), -rlim*cos(lon_rad), 'Color', [0 0 0], ...
    'LineWidth', 2);

axis tight
axis equal
hold(ax, 'off')
set(ax, 'FontSize', 12, 'Color', 'none', ...
    'XLim', [-1 1] * 1.02 * rlim, ...
    'YLim', [-1 1] * 1.02 * rlim)
set(get(ax, 'XAxis'), 'Visible', 'off')
set(get(ax, 'YAxis'), 'Visible', 'off')

varns = {ax};
varargout = varns(1:nargout);
end