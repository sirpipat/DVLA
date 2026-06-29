function varargout = realsot(sigma_n, sigma_m, lambda, ...
    sigma_trend, sigma_annual, sigma_semi)
% varargout = realsot(sigma_n, sigma_m, lambda, sigma_trend, ...
%     sigma_annual, sigma_semi)
%
% Calculates the time series data of the T-wave travel-time anomaly from 
% the correlated T-wave travel-time anomaly among repeaters (event pairs 
% that are close to each other and produce similar waveforms) which can be
% considered as the difference of the T-wave travel-time anomaly between
% the event pairs. It solves a linear equations by either assuming the
% linear, annual, and semi-annual trends or not.
%
% The inversion follows the method from Wu et al. (2023)
%
% INPUTS:
% sigma_n       prior uncertainty of the measurement nosise [default: 0.01]
% sigma_m       prior uncertainty of the model              [default: 0.15]
% lambda        smoothness parameter (in time-domain (days)). Set it to
%               zero to disable smoothing)                  [default: 15]
% sigma_trend   prior uncertainty of the coeff of the linear trend
%               [default: 15e-3 / 365.2425]
% sigma_annual  prior uncertainty of the coeff of the annual trend
%               [default: 15e-3]
% sigma_semi    prior uncertainty of the coeff of the semi-annual trend
%               [default: 15e-3]
%
% OUTPUT:
% tau_est       estimated T-wave travel-time anomaly for each event when no
%               the trend is assumed
% a_T_est       estimated stochastic part of the T-wave travel-time anomaly
%               for each event when the trend is assumed
% a_d           estimated coefficients of the deterministic part of the
%               T-wave travel-time anomaly for each event given by the
%               trend
% t             event origin time in days relative to the middle of the
%               time series
% dtm           datetime of the middle of the time series
%
% SEE ALSO:
% FAKESOT
%
% Last modified by spipatprathanporn@ucsd.edu, 06/29/2026
defval('sigma_n', 0.05)
defval('sigma_m', 0.15)
defval('lambda', 15)
defval('sigma_trend', 4e-3 / 356.2425)
defval('sigma_annual', 15e-3)
defval('sigma_semi', 15e-3)

omega = 2*pi/365.2425;

fname = fullfile('/Users/spipatprathanporn/research/IFILES/MATFILES', ...
    'sot_sumatra-psi-dgar_20260509.mat');
load(fname, 'cc_P', 'cc_T', 'dd', 'evs_str' , 'tt_shift_T');
A = and(and(cc_T >= 0.6, cc_P >= 0.9), dd <= 60);
wh = any(A,1);
A = A(wh', wh);
tt_shift_T = tt_shift_T(wh', wh);
m = sum(sum(A))/2;
n = size(A,1);

dt = datetime(evs_str.PreferredTime(wh)');
t = days(dt - datetime(2009,0,0));
tm = (t(end) + t(1)) / 2;
dtm = datetime(2009,0,0) + days(tm);
t = t - tm;


% construct a graph
g = graph(A);

% construct the linear equation
Eo = zeros(m,n);
delt = zeros(m,1);
for ii = 1:height(g.Edges)
    endnodes = g.Edges(ii,:).EndNodes;
    Eo(ii, endnodes(1)) = -1;
    Eo(ii, endnodes(2)) = 1;
    delt(ii) = -tt_shift_T(endnodes(1), endnodes(2));
end

%% Inversion
% covariance matrix of data
N = sigma_n^2 * eye(length(delt));

% prior covariance matrix of model parameters
R = 7.8503*sigma_m^2*7.8503 * ones(length(t));
if lambda
    R = R .* exp(-abs(t-t')/lambda);
else
    R = R.* eye(size(R));
end

% initial guess
tau0 = 0*t + 0;

% inversion
tau_est = (inv(R) + Eo'/N*Eo) \ (Eo'/N*delt + R\tau0);

% inversion with coefficients for deterministic part of T-wave anomaly
D = [t sin(omega*t) cos(omega*t) sin(2*omega*t) cos(2*omega*t)];
E = [Eo Eo*D];
S = blkdiag(R, ...
    7.8503*sigma_trend^2 * 7.8503, ...
    7.8503*sigma_annual^2 * 7.8503 / 2 * eye(2), ...
    7.8503*sigma_semi^2 * 7.8503 / 2 * eye(2));
a_est = (inv(S) + E'/N*E) \ (E'/N*delt);
% stochastic T-wave travel time anomaly
a_T_est = a_est(1:end-5);
% deterministic T-wave trave time anomaly
a_d = a_est(end-4:end);


%% Plot, I guess
% location of the earthquakes
x = evs_str.PreferredLongitude(wh);
y = evs_str.PreferredLatitude(wh);

% continuous time for the trend
t_truth = (0:0.002:1) * (t(end)-t(1)) + t(1);

% plot
figure(3)
clf
set(gcf, 'Units', 'inches', 'Position', [10 1 8 9])

% plot map
subplot(3,2,1)
hold on
for ii = 1:m
    endnodes = g.Edges(ii,:).EndNodes;
    plot(x(endnodes), y(endnodes), 'Color', 'k', ...
        'LineWidth', 0.05)
end
scatter(x,y,6,'y','filled','MarkerEdgeColor','k')
grid on
axis equal
xlabel('longitude')
ylabel('latitude')
title(sprintf('Map: %d events, %d repeaters', n, m))
set(gca, 'FontSize', 10, 'TickDir', 'out', 'Box', 'on')

% plot adjacency graph
subplot(3,2,2)
plot(g, 'Layout', 'force')
set(gca, 'FontSize', 10, 'TickDir', 'out', 'Box', 'on')
title(sprintf('Graph: %d events, %d repeaters', n, m))
axis tight

subplot(6,1,3)
hold on
for ii = 1:m
    endnodes = g.Edges(ii,:).EndNodes;
    plot(dt(endnodes), tau_est(endnodes), 'Color', rgbcolor('2'), ...
        'LineWidth', 0.05)
end
scatter(dt, tau_est, 6, rgbcolor('2'), "filled", 'o')
grid on
nolabels(gca, 1)
%ylabel('travel time anomaly (s)')
title('best fit without assuming trend')
set(gca, 'FontSize', 10, 'TickDir', 'out', 'Box', 'on')

subplot(6,1,4)
plot(datetime(2009,0,0) + days(t_truth + tm), ...
    Twavetimedelay(t_truth, [a_d' 0]), 'LineWidth', 1, 'Color', 'k')
hold on
scatter(dt, a_T_est + Twavetimedelay(t, [a_d' 0]), 6, rgbcolor('1'), ...
    "filled", 'o')
grid on
nolabels(gca, 1)
%ylabel('travel time anomaly (s)')
title(sprintf(['best fit with trends: (%.2g)t+%.2gsin(\\omegat)+' ...
    '%.2gcos(\\omegat)+%.2gsin(2\\omegat)+%.2gcos(2\\omegat)+' ...
    'a^{(T)}_{est}'], a_d(1), a_d(2), a_d(3), a_d(4), a_d(5)))
legend('deterministic (trend)', 'best fit')
set(gca, 'FontSize', 10, 'TickDir', 'out', 'Box', 'on')

subplot(6,1,5)
scatter(dt, a_T_est + Twavetimedelay(t, [a_d' 0]) - tau_est, 6, ...
    rgbcolor('1'), "filled", 'o')
grid on
nolabels(gca, 1)
ylabel('travel time anomaly (s)')
title('best fit with trends MINUS best fit without trends')
set(gca, 'FontSize', 10, 'TickDir', 'out', 'Box', 'on')

subplot(6,1,6)
scatter(dt, a_T_est, 6, rgbcolor('1'), 'filled')
grid on
xlabel('origin time')
%ylabel('travel time anomaly (s)')
title('stochastic T-wave anomaly (from inversion): a^{(T)}_{est}')
set(gca, 'FontSize', 10, 'TickDir', 'out', 'Box', 'on')

set(gcf, 'Renderer', 'painters')

% Matrices
figure(4)
clf
set(gcf, 'Units', 'inches', 'Position', [1 1 8 10])
subplot('Position', [0.07 0.96 0.88 0.01])
title(['Visualization: \delta = E_oɑ+n, ɑ_{est} = (R^{-1}+' ...
    'E_o^TN^{-1}E_o)^{-1}E_o^TN^{-1}\delta, E_o^TE_o = D-A'])
set(gca, 'FontSize', 12, 'Color', 'none')
nolabels(gca, 4)
set(get(gca, 'XAxis'), 'Visible', 'off')
set(get(gca, 'YAxis'), 'Visible', 'off')

subplot('Position', [0.07 0.695 0.41 0.26])
imagesc(A)
colormap(gca, [0.21 0.12 0.60; 0.9 0.75 1])
cb = colorbar;
set(cb, 'Limit', [-0.5 1.5], 'Ticks', [0 1], ...
    'TickLabels', {'no', 'yes'}, 'TickDirection', 'none')
set(get(cb, 'Label'), 'String', 'Are they adjacent?', 'FontSize', 11)
grid on
axis tight
axis equal
ylabel('index of first event')
xlabel('index of second event')
title('A: adjacency matrix')
set(gca, 'FontSize', 10, 'TickDir', 'out', 'Box', 'on')

subplot('Position', [0.57 0.695 0.41 0.26])
imagesc(Eo)
colormap(gca, [0.75 0 0; 1 1 0.9; 0 0 0.75])
clim([-1.5 1.5])
cb = colorbar ;
set(cb, 'Limit', [-1.5 1.5], 'Ticks', -1:1, 'TickDirection', 'none')
grid on
axis tight
axis equal
ylabel('index of repeater')
xlabel('index of event')
title('E_o: incidence matrix')
set(gca, 'FontSize', 10, 'TickDir', 'out', 'Box', 'on')

subplot('Position', [0.07 0.37 0.41 0.26])
imagesc(R)
colorbar
grid on
axis tight
axis equal
ylabel('index of first event')
xlabel('index of second event')
title('R: prior covariance of model')
set(gca, 'FontSize', 10, 'TickDir', 'out', 'Box', 'on')

subplot('Position', [0.57 0.37 0.41 0.26])
imagesc(N)
colorbar
grid on
axis tight
axis equal
ylabel('index of repeater')
xlabel('index of repeater')
title('N: prior covariance of data')
set(gca, 'FontSize', 10, 'TickDir', 'out', 'Box', 'on')

subplot('Position', [0.07 0.045 0.41 0.26])
imagesc(inv(R))
colorbar
grid on
axis tight
axis equal
ylabel('index of first event')
xlabel('index of second event')
title('R^{-1}')
set(gca, 'FontSize', 10, 'TickDir', 'out', 'Box', 'on')

subplot('Position', [0.57 0.045 0.41 0.26])
imagesc(Eo'/N*Eo)
colorbar
cm = colormap(gca);
cm(1,:) = [0.9 0.75 1];
colormap(gca, cm)
grid on
axis tight
axis equal
ylabel('index of event')
xlabel('index of event')
title('E_o^TN^{-1}E_o')
set(gca, 'FontSize', 10, 'TickDir', 'out', 'Box', 'on')

set(gcf, 'Renderer', 'painters')

% Matrices: Full version
figure(5)
clf
set(gcf, 'Units', 'inches', 'Position', [1 1 8 10])
subplot('Position', [0.07 0.96 0.88 0.01])
title(['Visualization: \delta = Eɑ+n, ɑ_{est} = (S^{-1}+' ...
    'E^TN^{-1}E)^{-1}E^TN^{-1}\delta'])
set(gca, 'FontSize', 12, 'Color', 'none')
nolabels(gca, 4)
set(get(gca, 'XAxis'), 'Visible', 'off')
set(get(gca, 'YAxis'), 'Visible', 'off')

subplot('Position', [0.07 0.695 0.41 0.26])
imagesc(A)
colormap(gca, [0.21 0.12 0.60; 0.9 0.75 1])
cb = colorbar;
set(cb, 'Limit', [-0.5 1.5], 'Ticks', [0 1], ...
    'TickLabels', {'no', 'yes'}, 'TickDirection', 'none')
set(get(cb, 'Label'), 'String', 'Are they adjacent?', 'FontSize', 11)
grid on
axis tight
axis equal
ylabel('index of first event')
xlabel('index of second event')
title('A: adjacency matrix')
set(gca, 'FontSize', 10, 'TickDir', 'out', 'Box', 'on')

subplot('Position', [0.57 0.695 0.41 0.26])
imagesc(E)
%clim([-1.5 1.5])
cb = colorbar ;
set(cb, 'TickDirection', 'out')
grid on
axis tight
axis equal
ylabel('index of repeater')
xlabel('index of event')
title('E = E_o [I D]')
set(gca, 'FontSize', 10, 'TickDir', 'out', 'Box', 'on')

subplot('Position', [0.07 0.37 0.41 0.26])
imagesc(S)
colorbar
grid on
axis tight
axis equal
ylabel('index of first event')
xlabel('index of second event')
title('S: prior covariance of model')
set(gca, 'FontSize', 10, 'TickDir', 'out', 'Box', 'on')

subplot('Position', [0.57 0.37 0.41 0.26])
imagesc(N)
colorbar
grid on
axis tight
axis equal
ylabel('index of repeater')
xlabel('index of repeater')
title('N: prior covariance of data')
set(gca, 'FontSize', 10, 'TickDir', 'out', 'Box', 'on')

subplot('Position', [0.07 0.045 0.41 0.26])
imagesc(inv(S))
colorbar
grid on
axis tight
axis equal
ylabel('index of first event')
xlabel('index of second event')
title('S^{-1}')
set(gca, 'FontSize', 10, 'TickDir', 'out', 'Box', 'on')

subplot('Position', [0.57 0.045 0.41 0.26])
imagesc(E'/N*E)
colorbar
grid on
axis tight
axis equal
ylabel('index of event')
xlabel('index of event')
title('E^TN^{-1}E')
set(gca, 'FontSize', 10, 'TickDir', 'out', 'Box', 'on')

set(gcf, 'Renderer', 'painters')

outputs = {tau_est, a_T_est, a_d, t, dtm};
varargout = outputs(1:nargout);
end

% define the absolute function here
function tau = Twavetimedelay(t,a)
    defval('a', [-1e-4 1/5 1/3 -1/8 -1/5 0])

    % period
    T = 365.2425;
    omega = 2*pi/T;

    tau = a(1)*t + a(2)*sin(omega*t) + a(3)*cos(omega*t) + ...
        a(4)*sin(2*omega*t) + a(5)*cos(2*omega*t) + a(6)*randn(size(t));
end