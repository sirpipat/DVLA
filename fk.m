function [ZZ, F, K] = fk(zz, dt, dx)
% ZZ = FK(zz, dt, dx)
%
% Perform f-k transformation on equally-spaced time-series data from an
% array on a straight line with an equal spacing.
%
% INPUT:
% zz        2D matrix representing time-series data from an array: 
%           zz = zz(t, x)
% dt        time spacing
% dx        spatial spacing
%
% OUTPUT:
% ZZ        2D Fourier transform of zz: ZZ = ZZ(F, K)
% F         frequency vector
% K         wave-number vector
%
%
% EXAMPLE:
% % this script
% t = (-10:0.01:20)'; dt = t(2) - t(1);
% x = 0:10:1000; dx = x(2) - x(1);
% [xx, tt] = meshgrid(x, t);
% zz = exp(-((tt - xx/250) / 0.25).^2) + 0.2 * rand(size(xx));
% [XX, F, K] = fk(xx, dt, dx);
%
% % generate this result
% fk('demo1');
%
% Last modified by spipatprathanporn@ucsd.edu, 11/24/2025

% Demo 1: simple wavefront
if strcmp(zz, 'demo1')
    t = (-10:0.01:20)'; dt = t(2) - t(1);
    x = 0:10:1000; dx = x(2) - x(1);
    [xx, tt] = meshgrid(x, t);
    zz = exp(-((tt - xx/250) / 0.25).^2) + 0.2 * rand(size(xx));
    [ZZ, F, K] = fk(zz, dt, dx);
    figure
    set(gcf, 'Units', 'inches', 'Position', [0 1 10 5])
    ax1 = subplot('Position', [0.07 0.12 0.38 0.80]);
    hold on
    for ii = 1:length(x)
        plot(t, zz(:,ii) * dx + x(ii), 'Color', 'k', 'LineWidth', 1);
    end
    grid on
    xlabel('time (s)')
    ylabel('location (m)')
    title('time-series plot (c = 250 m/s)')
    ylim([-10 1050])
    set(ax1, 'Box', 'on', 'TickDir', 'out')
    
    ax2 = subplot('Position', [0.53 0.12 0.44 0.80]);
    imagesc(K, F, 10*log10(abs(ZZ).^2)); 
    cb = colorbar; 
    grid on
    cb.Label.String = '10 log_{10}(spectral density)';
    xlabel('wave number (m^{-1})')
    ylabel('frequency (Hz)')
    title('F-K plot')
    set(ax2, 'Box', 'on', 'TickDir', 'out')

% Demo 2: several wavefronts with some dispersion
elseif strcmp(zz, 'demo2')
    t = (-10:0.01:20)'; dt = t(2) - t(1);
    x = 0:10:1000; dx = x(2) - x(1);
    [xx, tt] = meshgrid(x, t);
    zz = exp(-((tt - xx/250) / 0.25).^2) + 0.2 * rand(size(xx));
    
    for jj = 0:1200
        c = 1000 -jj/5;
        f = 2 + jj * 0.01;
        p = jj * pi/700;
        a = 1 - (jj - 600)/1000;
        ze = 1./(1200-xx) .* exp(-(tt - (1000-xx)/c - 3) / 4) ./ ...
            (1 + exp(-(tt - (1000-xx)/c - 3) / 0.1));
        zw = cos(2*pi*(tt - (1000-xx)/c - 3) * f + p) * a;
        zz = zz + ze .* zw;
    end

    for jj = 0:1200
        c = 2000 -jj/5;
        f = 2 + jj * 0.01;
        p = jj * pi/700;
        a = 1 - (jj - 600)/1000;
        ze = 1.2./(1200-xx) .* exp(-(tt - (1000-xx)/c - 8) / 4) ./ ...
            (1 + exp(-(tt - (1000-xx)/c - 8) / 0.1));
        zw = cos(2*pi*(tt - (1000-xx)/c - 8) * f + p) * a;
        zz = zz + ze .* zw;
    end

    [ZZ, F, K] = fk(zz, dt, dx);
    figure
    set(gcf, 'Units', 'inches', 'Position', [0 1 10 5])
    ax1 = subplot('Position', [0.07 0.12 0.38 0.80]);
    hold on
    for ii = 1:length(x)
        plot(t, zz(:,ii) * dx + x(ii), 'Color', 'k', 'LineWidth', 1);
    end
    grid on
    xlabel('time (s)')
    ylabel('location (m)')
    title('time-series plot')
    ylim([-10 1050])
    set(ax1, 'Box', 'on', 'TickDir', 'out')
    
    ax2 = subplot('Position', [0.53 0.12 0.44 0.80]);
    imagesc(K, F, 10*log10(abs(ZZ).^2)); 
    cb = colorbar; 
    grid on
    cb.Label.String = '10 log_{10}(spectral density)';
    xlabel('wave number (m^{-1})')
    ylabel('frequency (Hz)')
    title('F-K plot')
    set(ax2, 'Box', 'on', 'TickDir', 'out')
end

ZZ = fftshift(fft2(zz));

% list of frequencies
L = size(zz, 1);
if mod(L,2) == 0
    F = 1/dt/L * (-L/2:L/2-1);
else
    F = 1/dt/L * (-(L-1)/2:(L-1)/2);
end

% list of wave numbers
M = size(zz, 2);
if mod(M,2) == 0
    K = 2*pi/dx/M * (-M/2:M/2-1);
else
    K = 2*pi/dx/M * (-(M-1)/2:(M-1)/2);
end
end