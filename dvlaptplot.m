function figs = dvlaptplot(z, t, data, c, tipe)
% fig = DVLAFKPLOT(z, t, data, c, tipe)
%
% Makes a horizontal slowness (p) - time (t) plot for a vertical hydrophone
% array.
%
% INPUT:
% z         array hydrophone depths
% t         time
% data      hydrophone data: data = data(z, t)
% c         wave speed in the medium
% tipe      either
%           1 - vertical slowness stacking [default]
%           2 - horiztonal slowness stacking
%           3 - horizontal arriving angle
%           4 - incident angle
%
% OUTPUT:
% figs       figure handles
%
% Last modified by spipatprathanporn@ucsd.edu, 04/30/2026

% tipe
% 1 - vertical slowness (dt/dz)
% 2 - horizontal slowness (sqrt(1/c^2-(dt/dz)^2))
% 3 - horizontal arriving angle (asin(c * dt/dz))
% 4 - incident angle (acos(c * dt/dz))
defval('tipe', 1)

if tipe == 1
    % vertical slowness
    eta = (-1:0.002:1)' * 1/c;
    % horizontal slownesses
    p = sqrt(1/c^2 - eta.^2);
    yval = eta;
    ytext = 'vertical slowness (s/m)';
elseif tipe == 2
    % horizontal slownesses
    p = (-1:0.002:1)' * 1/c;
    % vertical slowness
    eta = sqrt(1/c^2 - p.^2) .* (-1 + 2 * (p >= 0));
    yval = p;
    ytext = 'horizontal slowness (s/m)';
elseif tipe == 3
    % horizontal arriving angle
    theta = (-1:0.001:1)' * pi/2;
    % vertical slowness
    eta = sin(theta) / c;
    % horizontal arriving angle in degrees
    yval = theta * 180/pi;
    ytext = 'horizontal arriving angle (deg)';
elseif tipe == 4
    % incident angle
    theta = (-1:0.001:1)' * pi/2;
    % vertical slowness
    eta = cos(theta) / c;
    % horizontal arriving angle in degrees
    yval = theta * 180/pi;
    ytext = 'incident angle (deg)';
else
    error('tipe must be 1, 2, 3, or 4.')
end

% moving window is 0.1 seconds long
fs = 1 / (t(2) - t(1));
N = round(1*fs);

pt = zeros([length(eta) length(t)]);

% Stacked seismograms for plotting
% The number of stacked seismograms are reduced by a factor of 10 to keep
% the number managable.
%zs = zeros([length(yval) length(t)]);
zs = zeros(size(pt));

% Stack the seismograms
for jj = 1:length(eta)
    xs = data(1,:) / median(abs(data(1,:)));
    for ii = 2:length(z)
        xs = xs + interp1(t + eta(jj) * (z(ii) - z(1)), ...
            data(ii,:) / max(median(abs(data(ii,:))), 1), t, 'linear', 0);
    end
    xs_mvar = movvar(xs, N);
    pt(jj,:) = xs_mvar;
    zs(jj, :) = xs;
end
zs_scaling = max(abs(zs), [], 'all') / (yval(101)-yval(1)) * 0.8;

% autocorrelation of stacked seismograms
% azs = zeros([size(zs,1) size(zs,2)*2-1]);
% for ii = 1:length(eta)
%     azs(ii,:) = xcorr(detrend(pt(ii,:), 0), 'coeff');
% end

fig1 = figure(1);
clf
set(gcf, 'Units', 'inches', 'Position', [0 3 8 6])
hold on
for ii = 1:100:length(yval)
    plot(t, zs(ii,:) / zs_scaling + yval(ii), 'LineWidth', 1, 'Color', 'k')
end
grid on
set(gca, 'FontSize', 12, 'TickDir', 'out', 'Box', 'on')
axis tight
xlabel('time (s)')
ylabel(ytext)
set(gcf, 'Renderer', 'painters')

fig2 = figure(2);
clf
set(gcf, 'Units', 'inches', 'Position', [0 1 8 6])
imagesc([t(1) t(end)], [yval(1) yval(end)], pt.^0.5)
axis xy
grid on
set(gca, 'FontSize', 12, 'TickDir', 'out', 'Box', 'on')
xlabel('time (s)')
ylabel(ytext)

colormap('jet')
cb = colorbar;
set(cb.Label, 'String', '1-s moving std')
set(cb, 'TickDirection', 'both')
%set(gca, "ylim", [0 1] .* ylim)
set(gcf, 'Renderer', 'painters')



% fig4 = figure(4);
% clf
% set(gcf, 'Units', 'inches', 'Position', [0 5 8 6])
% imagesc([t(1) t(end)], [yval(1) yval(end)], abs(zs).^0.5)
% axis xy
% grid on
% set(gca, 'FontSize', 12, 'TickDir', 'out', 'Box', 'on')
% xlabel('time (s)')
% ylabel(ytext)
% 
% colormap('jet')
% cb = colorbar;
% set(cb.Label, 'String', 'sqrt(|pressure|)')
% set(cb, 'TickDirection', 'both')
% %set(gca, "ylim", [0 1] .* ylim)
% set(gcf, 'Renderer', 'painters')

% fig5 = figure(5);
% clf
% set(gcf, 'Units', 'inches', 'Position', [0 7 8 6])
% imagesc((t(end)-t(1)) * [-1 1], [yval(1) yval(end)], azs)
% axis xy
% grid on
% set(gca, 'FontSize', 12, 'TickDir', 'out', 'Box', 'on')
% xlabel('time (s)')
% ylabel(ytext)

% colormap('kelicol')
% clim([-1 1])
% %xlim([-10 10])
% cb = colorbar;
% set(cb.Label, 'String', 'autocorrelation of pressure')
% set(cb, 'TickDirection', 'both')
% set(gcf, 'Renderer', 'painters')


figs = [fig1 fig2];
end