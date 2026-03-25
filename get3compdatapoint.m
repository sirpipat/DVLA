function get3compdatapoint(src, ax, data, name, ddir, evid)
% GET3COMPDATAPOINT(src, ax, data, name, ddir, evid)
%
% Display a data point on the scatter plot. Call this function to the 
% scatter plot's button down function. Then, click one of the point to get
% the (x,y) coordinate as well as the index of the data point. If DATA is
% specified, the value at the index will be displayed as well. In addition,
% the seismogram plots from PLOT3COMPONENTS will show up as well if they
% alreay exist.
%
% INPUT:
% src       selected object
% ax        parent axes of the selected object
% data      data structure containing the plotted data [optional]
%           it could be a matrix, cell array, or struct.
% name      name of the data [optional]
% ddir      directory to the printed figures from PLOT3COMPONENTS
%
% SEE ALSO:
% GETDATAPOINT
%
% Last modified by spipatprathanporn@ucsd.edu, 03/26/2026

defval('data', [])
defval('name', 'value')

getdatapoint(src, ax, data, name)

pt = ax.CurrentPoint(1,1:2);
if ~(isempty(src.XData) || isempty(src.YData))
    ii = knnsearch([src.XData' src.YData'], pt);
else
    % Try using polar coordinates
    % Make sure that theta falls between 0 and 2 * pi radians
    ii = knnsearch([mod(src.ThetaData', 2*pi) src.RData'], pt);
end
if ~isempty(src.UserData)
    ii = src.UserData(ii);
end
fprintf('x : %g\n', pt(1));
fprintf('y : %g\n', pt(2));
fprintf('ii: %g\n', ii);

net = strip(data.KNETWK{ii});
sta = strip(data.KSTNM{ii});

system(sprintf('open -a Preview %s/plot3components_%d_%s.%s.*.pdf', ...
    ddir, evid, net, sta));
end