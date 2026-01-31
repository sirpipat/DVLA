function xcorrdvla(station, fs, fc, dt_begin, dt_end, scaling)
% XCORRDVLA(station, fs, fc, dt_begin, dt_end, scaling)
%
% Make the auto- and cross- correlation plots of the hydroacoustic
% seismograms recorded by a vertical hydrophone array.
%
% INPUT:
% station       station name
% fs            sampling rate
% fc            corner frequency for bandpass filtering
% dt_begin      beginning datetime
% dt_end        ending datetime
% scaling       amplitude scaling for plotting seismograms
%
% Last modified by spipatprathanporn@ucsd.edu, 01/30/2026

defval('fs', 10)
defval('fc', [0.5 5])
defval('scaling', 1)

% list the hydrophones
[list, depth] = getdvlaseis(station, 'list');
dz = depth(2) - depth(1);

% determine sampling times
dt = (dt_begin:seconds(1/fs):dt_end);

% construct the seismic data array
xx = nan(length(list), length(dt));

% read the seismogram
for ii = 1:length(list)
    try
        [t, x, jdn] = getdvlaseis(station, list(ii), ...
            dt_begin - seconds(1), dt_end + seconds(1));
        dt_dvla = datetime(jdn, "ConvertFrom", "datenum", ...
            "Format", "uuuu-MM-dd'T'HH:mm:ss.SSS") + seconds(t);
        fs0 = (length(t) - 1) / (t(end) - t(1));
    
        % lowpass filter to remove the potential alias
        x = lowpass(x, fs0, fs/2, 2, 1, 'butter', 'linear');

        % resample to the target sampling rate
        x = interp1(dt_dvla, x, dt);

        % bandpass filtering
        xx(ii,:) = bandpass(x, fs, fc(1), fc(2), 2, 1, 'butter', 'linear');
    catch ME
        ME.getReport()
        continue
    end
end

% plot seismograms
figure(10)
set(gcf, "Units", "inches", "Position", [0 1 10 8]);
clf
ax1 = subplot("Position", [0.07 0.57 0.40 0.40]);
hold on
x_scaling = max(abs(xx), [], "all");
for ii = 1:length(list)
    plot(dt, -xx(ii,:) / x_scaling * dz * scaling + depth(ii), ...
        'LineWidth', 1, 'Color', 'k')
end
axis ij
axis tight
xlabel("time")
ylabel("hydrophone depth (m)")
title(sprintf("%s (%d): %g-%g Hz", station, ...
    floor(jdn - datenum(2018, 12, 31)), fc(1), fc(2)))
grid on
set(ax1, "FontSize", 12, "Box", "on", "TickDir", "out")


% plot autocorrelations
ax2 = subplot("Position", [0.57 0.57 0.40 0.40]);
hold on
for ii = 1:length(list)
    % compute autocorrelation
    [r, lag] = xcorr(xx(ii,:), "coeff");
    % Note: the negative sign corrects the sign when plotting with depth 
    plot(lag/fs, -r * dz/1.5 + depth(ii), 'LineWidth', 1, 'Color', 'k')
end
axis ij
axis tight
xlabel("lag (s)")
ylabel("hydrophone depth (m)")
title("autocorrelation")
grid on
xlim([-10 10])
set(ax2, "FontSize", 12, "Box", "on", "TickDir", "out")


% plot correlation with the top most hydrophone
rr = nan(length(list), 2*length(dt)-1);

ax3 = subplot("Position", [0.07 0.07 0.40 0.40]);
hold on
for ii = 1:length(list)
    % compute correlation with the top most hydrophone
    [r, lag] = xcorr(xx(1,:), xx(ii,:), "coeff");
    rr(ii,:) = r;
    % Note: the negative sign corrects the sign when plotting with depth 
    plot(lag/fs, -r * dz/1.5 + depth(ii), 'LineWidth', 1, 'Color', 'k')
end
axis ij
axis tight
xlabel("lag (s)")
ylabel("hydrophone depth (m)")
title("correlation with the top most hydrophone")
grid on
xlim([-10 10])
set(ax3, "FontSize", 12, "Box", "on", "TickDir", "out")

% stack
% only consider the positive timeshift
% TODO: FK Filtering to remove sea-surface reflection
rr_pos = rr(:,length(dt):end);
lag_pos = lag(length(dt):end);
[MM, II] = max(rr_pos, [], 2);
timelag = lag_pos(II) / fs;

% determine slowness
% p(1) = slowness, p(2) = y-intercept of timelag-depth plot
p = polyfit(depth - depth(1), timelag, 1);
p = abs(p);

% plot what is going to be stacked
ax4 = subplot("Position", [0.57 0.07 0.40 0.40]);
hold on
x_scaling = max(abs(xx), [], "all");
for ii = 1:length(depth)
    plot(dt + seconds(p(1) * (depth(ii) - depth(1))), ...
        -xx(ii,:) / x_scaling * dz * scaling + depth(ii), ...
        'LineWidth', 1, 'Color', 'k')
end
axis ij
axis tight
xlabel("time")
ylabel("hydrophone depth (m)")
title(sprintf("Slowness stacking: slowness = %.3e s/m", p(1)))
grid on
set(ax4, "FontSize", 12, "Box", "on", "TickDir", "out")
set(gcf, "Renderer", "painters")

% plot stacked waveform
x_stack = xx(1,:);
for ii = 2:length(list)
    x_stack = x_stack + interp1(dt + seconds(p(1) * ...
        (depth(ii) - depth(1))), xx(ii,:), dt, 'linear', 0);
end
plot(dt, -x_stack / abs(max(x_stack)) * dz + (depth(1) - 2*dz), ...
        'LineWidth', 1.5, 'Color', 'r')
xlim([dt(1) dt(end)])
end