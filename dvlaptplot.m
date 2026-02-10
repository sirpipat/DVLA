function fig = dvlaptplot(z, t, data, c)
% fig = DVLAFKPLOT(z, t, data, c)
%
% Makes a horizontal slowness (p) - time (t) plot for a vertical hydrophone
% array.
%
% INPUT:
% z         array hydrophone depths
% t         time
% data      hydrophone data: data = data(z, t)
% c         wave speed in the medium
%
% OUTPUT:
% fig       figure handle
%
% Last modified by spipatprathanporn@ucsd.edu, 02/09/2026

% horizontal slownesses
p = (-1:0.01:1)' * 1/c;
% vertical slowness
eta = sqrt(1/c^2 - p.^2);

% moving window is 10 seconds long
fs = 1 / (t(2) - t(1));
N = round(10*fs);

pt = zeros([length(p) length(t)]);

for jj = 1:length(p)
    xs = data(1,:);
    for ii = 2:length(z)
        xs = xs + interp1(t + sign(p(jj)) * ...
            eta(jj) * (z(ii) - z(1)), ...
            data(ii,:), t, 'linear', 0);
    end
    xs_mvar = movvar(xs, N);
    pt(jj,:) = xs_mvar;
end

fig = figure(2);
clf
set(gcf, 'Units', 'inches', 'Position', [0 1 8 6])
imagesc([t(1) t(end)], [p(1) p(end)], log10(pt))
axis xy
grid on
set(gca, 'FontSize', 12, 'TickDir', 'out', 'Box', 'on')
xlabel('time (s)')
ylabel('horizontal slowness (s/m)')

cb = colorbar;
set(cb.Label, 'String', 'log_{10} 10-s moving variance')
set(cb, 'TickDirection', 'both')
set(gcf, 'Renderer', 'painters')
end