function varargout = plotpolarmap2(func_handl, ax, varargin)
% PLOTPOLARMAP2(lon, lat)
% PLOTPOLARMAP2(ax, lon, lat)
% PLOTPOLARMAP2(func_handl, lon, lat)
% PLOTPOLARMAP2(func_handl, ax, lon, lat)
% PLOTPOLARMAP2(..., varargin)
% obj = PLOTPOLARMAP2(...)
%
% Wrapper function for plotting anything on a polarmap made by
% POLARBASEMAP2.
%
% INPUT:
% func_handl    function handle for plotting    [default: @plot]
% ax            target axes                     [default: gca]
% lon           lontitudes
% lat           latitudes
% varargin      other input arguments for func_handl
%
% OUTPUT:
% obj           object handle
%
% EXAMPLE:
% ax = polarbasemap(40);
% sc = plotpolarmap(@scatter, [0 30 60 90], [50 60 70 80], 50, 'y', ...
%     'filled', 'v', 'MarkerEdgeColor', 'k');
%
% SEE ALSO:
% POLARBASEMAP2, PLOTPOLARMAP
%
% Last modified by spipatprathanporn@ucsd.edu, 08/10/2026

if isa(func_handl, 'function_handle')
    if ~isa(ax, 'matlab.graphics.axis.Axes')
        varargin = [{ax} varargin];
        ax = gca;
    end
else
    if ~isa(func_handl, 'matlab.graphics.axis.Axes')
        varargin = [{func_handl} {ax} varargin];
        func_handl = @plot;
        ax = gca;
    else
        varargin = [{ax} varargin];
        ax = func_handl;
        func_handl = @plot;
    end
end

% compute (x,y) coordinates from (lon,lat)
lon = varargin{1};
lat = varargin{2};
[r, az] = distance(ax.UserData(2), ax.UserData(1), lat, lon, 'degrees');

x = r .* sin(az*pi/180);
y = r .* cos(az*pi/180);

% making plot
hold(ax, 'on')
obj = feval(func_handl, ax, x, y, varargin{3:end});

% gather the output
varns = {obj};
varargout = varns(1:nargout);
end