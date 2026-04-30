function ax = plotdvlaseisarray(station, dt_begin, dt_end, colo, cohi)
% ax = PLOTDVLASEISARRAY(station, dt_begin, dt_end, colo, cohi)
%
% Plots a short section of the records from a hydrophone mooring.
%
% INPUT:
% station       station name either 'NERSC1', 'NERSC2', 'NERSC3', 'SIO1', 
%               'SIO2', or 'SIO3'
% dt_begin      beginning datetime or datestr
% dt_end        ending datetime or datestr
% colo          lower corner frequency for filtering
% cohi          higher corner frequency for filtering
%
% OUTPUT:
% ax            axes handle of the plot
%
% Last modified by spipatprathanporn@ucsd.edu, 04/30/2026

figure(3);
set(gcf, "Units", "inches", "Position", [0 1 8 6]);
clf
hold on

% getting the list of the hydrophone number and depth
[list, depth] = getdvlaseis(station, 'list');

% hydrophone separation (for amplitude scaling in the plot)
dz = depth(2) - depth(1);

% making plot
for ii = 1:length(list)
    try
        % obtain a section of seismogram
        [t, x, jdn] = getdvlaseis(station, list(ii), dt_begin, dt_end);

        % downsampling to by a factor of 20
        fs = (length(t) - 1) / (t(end) - t(1));
        xf = lowpass(x, fs, 10, 2, 1, 'butter', 'linear');
        xd = xf(1:20:end);
        td = t(1:20:end);
        fsd = (length(td) - 1) / (td(end) - td(1));
        % bandpass filtering
        xdf = bandpass(xd, fsd, colo, cohi, 2, 1, 'butter', 'linear');

        % plot
        plot(td, xdf / max(abs(xdf)) * dz + depth(ii), 'Color', 'k', ...
            'LineWidth', 1);
    catch ME
        ME.getReport()
        continue
    end
end

axis tight
axis ij
grid on
box on
ax = gca;
set(ax, "TickDir", "out", "FontSize", 12)
xlabel(sprintf("time since %s (s)", datetime(jdn, "ConvertFrom", ...
    "datenum", "Format", "uuuu-MM-dd HH:mm:ss.SSS")));
ylabel("hydrophone depth (m)")
title(sprintf("%s (%d): %g-%g Hz", station, ...
    floor(jdn - datenum(2018, 12, 31)), colo, cohi))
set(gcf, 'Renderer', 'painters')
end