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
% Last modified by spipatprathanporn@ucsd.edu, 01/30/2026

% slownesses
p = (-1:0.01:1)' * 1/c;

pt = zeros([length(p) length(t)]);

for jj = 1:length(p)
    xs = data(1,:);
    for ii = 2:length(list)
        xs = xs + interp1(dt + seconds(sign(p(jj)) * ...
            sqrt(1/c^2 - pp(jj)^2) * (z(ii) - z(1))), ...
            data(ii,:), dt, 'linear', 0);
    end
    xs_mvar = movvar(xs, N);
    pt(jj,:) = xs_mvar;
end

fig = figure(2);
clf
set(gcf, 'Units', 'inches', 'Position', [0 1 8 6])
imagesc([t(1) t(end)], [p(1) p(end)], pt)
axis xy
grid on
set(gca, 'FontSize', 12, 'TickDir', 'out', 'Box', 'on')
xlabel('time (s)')
ylabel('horizontal slowness (s/m)')

cb = colorbar;
set(cb.Label, 'String', '10-s moving variance')
set(cb, 'TickDirection', 'both')
set(gcf, 'Renderer', 'painters')
end