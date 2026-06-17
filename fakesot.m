function varargout = fakesot(op1, op2, op3, a, sigma_n, sigma_m, lambda)
% varargout = fakesot(op1, op2, op3, [a1 a2 a3 a4 a5 a6], sigma_n, ...
%     sigma_m, lambda)
%
% Generates a fake sesimic ocean thermography experiment data set to test
% inversion method. It attempts to solve the simplified version of equation
% (2) from Wu et al. (2023).
%
% INPUT:
% op1       options for event locations and repeater connectivity
%   - 1       random event location, random repeater connectivity
%   - 2       random event location, "snowflakes" repeater connectivity
%   - other   sumatra events location, and repeater connectivity  [default]
% op2       options for repeater spreading in time
%   - 1       repeater tend to be well-sampled in time            [default]
%   - 2       repeater tend to cluster together with some mixing
%   - other   repeater subgraph are sequential (bad sampling in time)
% op3       options for origin time sampling
%   - 1       uniform sampling
%   - other   - If op1 == 1 or 2, the Omori-Utsu's law            [default]
%             - Otherwise       , sumatra event origin times      [default]
%
% [a1 a2 a3 a4 a5 a6]         model parameters for ground truth
%       function of the T-wave travel-time anomaly (tau) for each event 
%       tau = [t sin(w*t) cos(w*t) sin(2*w*t) cos(2*w*t) randn()] * a'
%       where w = 2*pi/(1 year). Please note that a1-a5 are for the 
%       deterministic portion while a6 is for stochaistic contribution of 
%       the T-wave travel-time anomaly for each event.
%       [default: [-1e-4 1/5 1/3 -1/8 -1/5 0]]
% sigma_n       prior uncertainty of the measurement nosise [default: 0.05]
% sigma_m       prior uncertainty of the model              [default: 1]
% lambda        smoothness parameter (in time-domain (days)). Set it to
%               zero to disable smoothing)                  [default: 15]
%
% OUTPUT:
% tau_est       estimated T-wave travel-time anomaly for individual events
% tau           ground truth T-wave travelt-time anomaly
% Eo            design matrix of the linear equation
% m             number of repeaters
% n             number of unique events
% N             covariance matrix of the data noise
% R             prior covariance matrix of the model (tau)
% A             adjacency matrix of a graph of events and repeaters
%
% Last modified by spipatprathanporn@ucsd.edu, 06/17/2026

defval('op1', 0)
defval('op2', 1)
defval('op3', 0)
defval('a', [-1e-4 1/5 1/3 -1/8 -1/5 0])
defval('sigma_n', 0.05)
defval('sigma_m', 1)
defval('lambda', 15)

%% Part 1: Generating a fake data set
T = 365.2425;
omega = 2*pi/T;
t_truth = (-3:0.01:3) * T;
tau_truth = Twavetimedelay(t_truth, a);

if op1 == 1
    % make some fictional data and model
    n = 0;
    while n <= 0
        % number of events
        n = 800;
        
        % begin and end of the time series
        t0 = -3*T;
        t1 =  3*T;
        
        % event origin times
        if op3 == 1
            t = sort(t0 + rand(n,1) * (t1 - t0));
        else
            t = sort(randeventtime(t0, t1, 2*T, 1, n, []));
        end
    
        % event hypocenter
        x0 = rand(n/8,3);
        x = x0 * [300/sqrt(2) -300/sqrt(2) 0; ...
                  100/sqrt(2)*cos(pi/18) 100/sqrt(2)*cos(pi/18) 100*sin(pi/18); ...
                  -30/sqrt(2)*sin(pi/18) -30/sqrt(2)*sin(pi/18) 30*cos(pi/18)];
        x = kron(x, ones(8,1)) + rand(n,3) * diag([30 30 6]);
    
        % sort the time
        % [t, it] = sort(t);
        % x = x(it,:);
    
        % event separations
        dd = sqrt((x(:,1)-x(:,1)').^2 + (x(:,2)-x(:,2)').^2 + (x(:,3)-x(:,3)').^2);
        
        % (absolute) T-wave travel time anomaly
        tau = Twavetimedelay(t, a);
        
        % cross-correlation table
        CC = rand(n,n);
        % make it symmetric
        CC = ((CC + CC') / 2);
        % make cc value stronger when the event separation is small
        CC = (CC.^1) .* exp(-dd/36);
        % make its diagonal zero since we are looking for event pairs
        CC = CC .* (1 - eye(size(CC)));
        
        % connectivity matrix
        A = and(CC >= 0.6, dd<=60);
        
        % number of repeaters
        m = sum(sum(A))/2;
        
        % remove unused events
        used_events = any(A,2);
        n = sum(used_events);
    end
    t = t(used_events);
    tau = tau(used_events);
    A = A(used_events, used_events');
    CC = CC(used_events, used_events');
    x = x(used_events,:);

    if op2 ~= 1
        [A, p] = forceblockdiagonal(A);
        CC = CC(p,p);
        x = x(p,:);
        if op2 == 2
            nperm = 20;
            ngroup = floor(size(A,1) / nperm);
            nremainder = mod(size(A,1), nperm);
            [~, idx] = sort(rand(nperm, ngroup), 1);
            iperm = reshape(idx + ones(nperm,1) * (0:ngroup-1) * nperm, 1, []);
            iperm = [iperm randperm(nremainder, nremainder)+ngroup*nperm];
            A = A(iperm, iperm);
            CC = CC(iperm, iperm);
            x = x(iperm,:);
        end
    end
elseif op1 == 2
    % number of events
    n = 320;
    
    % begin and end of the time series
    t0 = -3*T;
    t1 =  3*T;
    
    % event origin times
    if op3 == 1
        t = t0 + rand(n,1) * (t1 - t0);
    else
        t = randeventtime(t0, t1, 2*T, 1, n, []);
    end

    % event hypocenter
    x0 = rand(n/8,3);
    x = x0 * [300/sqrt(2) -300/sqrt(2) 0; ...
              100/sqrt(2)*cos(pi/18) 100/sqrt(2)*cos(pi/18) 100*sin(pi/18); ...
              -30/sqrt(2)*sin(pi/18) -30/sqrt(2)*sin(pi/18) 30*cos(pi/18)];
    x = kron(x, ones(8,1)) + rand(n,3) * diag([15 15 3]);

    % sort the time
    [t, it] = sort(t);
    x = x(it,:);

    % event separations
    dd = sqrt((x(:,1)-x(:,1)').^2 + (x(:,2)-x(:,2)').^2 + (x(:,3)-x(:,3)').^2);
    
    % (absolute) T-wave travel time anomaly
    tau = Twavetimedelay(t, a);

    % force the connectivity graph to take snowflakes shape
    [~,A0] = fakesnowflakes(n, floor(n/5.4), 2);
    if op2 == 1
        iperm = randperm(size(A0,1), size(A0,1));
        A0 = A0(iperm,iperm);
    elseif op2 == 2
        nperm = 20;
        ngroup = floor(size(A0,1) / nperm);
        nremainder = mod(size(A0,1), nperm);
        [~, idx] = sort(rand(nperm, ngroup), 1);
        iperm = reshape(idx + ones(nperm,1) * (0:ngroup-1) * nperm, 1, []);
        iperm = [iperm randperm(nremainder, nremainder)+ngroup*nperm];
        A0 = A0(iperm, iperm);
    end
    CC = A0 .* (0.5 + 0.5*randsym(size(A0))) + ~A0 .* randsym(size(A0)) .* ...
        (0.6-eps) .* ~eye(size(A0));

    % connectivity matrix
    A = and(CC >= 0.6, dd <= 60);

    % number of repeaters
    m = sum(sum(A))/2;

    % remove unused events
    used_events = any(A,2);
    n = sum(used_events);
    t = t(used_events);
    tau = tau(used_events);
    A = A(used_events, used_events');
    CC = CC(used_events, used_events');
    x = x(used_events,:);
else
    n = 203;
    % begin and end of the time series
    t0 = -3*T;
    t1 =  3*T;

    fname = fullfile('/Users/spipatprathanporn/research/IFILES/MATFILES', 'sot_sumatra-psi-dgar_20260509.mat');
    load(fname, 'cc_P', 'cc_T', 'dd', 'evs_str');
    A = and(and(cc_T >= 0.6, cc_P >= 0.9), dd <= 60);
    wh = any(A,1);
    A = A(wh', wh);
    m = sum(sum(A))/2;
    x = deg2km(evs_str.PreferredLongitude(wh)' - 95);
    y = deg2km(evs_str.PreferredLatitude(wh)' - 0.5);
    z = evs_str.PreferredDepth(wh)';
    x = [x y z];

    if op2 ~= 1
        [A, p] = forceblockdiagonal(A);
        x = x(p,:);
        if op2 == 2
            nperm = 20;
            ngroup = floor(size(A,1) / nperm);
            nremainder = mod(size(A,1), nperm);
            [~, idx] = sort(rand(nperm, ngroup), 1);
            iperm = reshape(idx + ones(nperm,1) * (0:ngroup-1) * nperm, 1, []);
            iperm = [iperm randperm(nremainder, nremainder)+ngroup*nperm];
            A = A(iperm, iperm);
            x = x(iperm,:);
        end
    end

    % event origin times
    if op3 == 1
        t = t0 + rand(n,1) * (t1 - t0);
        t = sort(t);
    else
        t = days(datetime(evs_str.PreferredTime(wh)') - datetime(2009,0,0));
    end
    % (absolute) T-wave travel time anomaly
    tau = Twavetimedelay(t, a);
    
end

% construct a graph
g = graph(A);

% construct the linear equation
Eo = zeros(m,n);
for ii = 1:height(g.Edges)
    endnodes = g.Edges(ii,:).EndNodes;
    Eo(ii, endnodes(1)) = -1;
    Eo(ii, endnodes(2)) = 1;
end

% observed data
delt = Eo * tau + 0. * randn(m,1);

%% Part 2: Inversion
% covariance matrix of data
N = sigma_n^2 * eye(length(delt));

% prior covariance matrix of model parameters
R = sigma_m^2 * ones(length(tau));
if lambda
    R = R .* exp(-abs(t-t')/lambda);
else
    R = R.* eye(size(R));
end

% initial guess
tau0 = 0*tau + 0;

% inversion
tau_est = (inv(R) + Eo'/N*Eo) \ (Eo'/N*delt + R\tau0);

% inversion with coefficients for deterministic part of T-wave anomaly
D = [t sin(omega*t) cos(omega*t) sin(2*omega*t) cos(2*omega*t)];
E = [Eo Eo*D];
S = blkdiag(R, ...
    7.8503*(4e-3 / 356.2425)^2 * 7.8503, ...
    7.8503*(15e-3)^2 * 7.8503 / 2 * eye(2), ...
    7.8503*(15e-3)^2 * 7.8503 / 2 * eye(2));
a_est = (inv(S) + E'/N*E) \ (E'/N*delt);
% stochastic T-wave travel time anomaly
a_T_est = a_est(1:end-5);
% deterministic T-wave trave time anomaly
a_d = a_est(end-4:end);

%% Part 3: Plot
% plot
figure(3)
clf
set(gcf, 'Units', 'inches', 'Position', [10 1 8 9])

% plot map
subplot(3,2,1)
hold on
for ii = 1:m
    endnodes = g.Edges(ii,:).EndNodes;
    plot(x(endnodes,1), x(endnodes,2), 'Color', 'k', ...
        'LineWidth', 0.05)
end
scatter(x(:,1),x(:,2),6,'y','filled','MarkerEdgeColor','k')
grid on
axis equal
%xlim([-50 350])
%ylim([-200 100])
xlabel('Easting (km)')
ylabel('Northing (km)')
title(sprintf('Map: %d events, %d repeaters', n, m))
set(gca, 'FontSize', 10, 'TickDir', 'out', 'Box', 'on')

% plot adjacency graph
subplot(3,2,2)
plot(g, 'Layout', 'force')
set(gca, 'FontSize', 10, 'TickDir', 'out', 'Box', 'on')
title(sprintf('Graph: %d events, %d repeaters', n, m))
axis tight

% plot ground truth repeaters
subplot(3,2,3)
plot(t_truth, tau_truth, 'LineWidth', 0.5, 'Color', [0.6 0.6 0.6])
hold on
for ii = 1:m
    endnodes = g.Edges(ii,:).EndNodes;
    plot(t(endnodes), interp1(t_truth, tau_truth, t(endnodes)), '-o', ...
        'Color', rgbcolor('1'), 'MarkerFaceColor', rgbcolor('1'), ...
        'MarkerSize', 2, 'LineWidth', 0.05)
end
grid on
xlabel('origin time (days)')
ylabel('travel time anomaly (s)')
title('repeaters in time series')
set(gca, 'FontSize', 10, 'TickDir', 'out', 'Box', 'on')
ylim([-1 1]*1)

subplot(3,2,4)
plot(t_truth, tau_truth, 'LineWidth', 1)
hold on
scatter(t, tau_est, 6, rgbcolor('2'), "filled", 'o')
grid on
xlabel('origin time (days)')
ylabel('travel time anomaly (s)')
title('best fit')
legend('ground truth', 'best-fit')
set(gca, 'FontSize', 10, 'TickDir', 'out', 'Box', 'on')
ylim([-1 1]*1+[0 0.5])

subplot(3,2,5)
plot(t, a_T_est, 'Color', rgbcolor('2'), 'LineWidth', 1)
grid on
xlabel('origin time (days)')
ylabel('travel time anomaly (s)')
title('stochastic T-wave anomaly (from inversion)')
set(gca, 'FontSize', 10, 'TickDir', 'out', 'Box', 'on')
if a(6) > 0
    ylim([-3 3]*abs(a(6)))
else
    ylim([-0.5 0.5])
end

subplot(3,2,6)
plot(t_truth, tau_truth, 'LineWidth', 1)
hold on
plot(t_truth, Twavetimedelay(t_truth, [a_d' 0]), 'LineWidth', 0.5)
scatter(t, a_T_est + Twavetimedelay(t, [a_d' 0]), 6, rgbcolor('2'), "filled", 'o')
grid on
xlabel('origin time (days)')
ylabel('travel time anomaly (s)')
title('best fit with trends')
legend('ground truth', 'deterministic (trend)', 'best fit')
set(gca, 'FontSize', 10, 'TickDir', 'out', 'Box', 'on')
ylim([-1 1]*1+[0 0.5])
set(gcf, 'Renderer', 'painters')

% Matrices
figure(4)
clf
set(gcf, 'Units', 'inches', 'Position', [1 1 8 10])
subplot('Position', [0.07 0.96 0.88 0.01])
title('Visualization: \delta = E_oɑ+n, ɑ_{est} = (R^{-1}+E_o^TN^{-1}E_o)^{-1}E_o^TN^{-1}\delta, E_o^TE_o = D-A')
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

outputs = {tau_est, tau, Eo, m, n, N, R, A};
varargout = outputs{1:nargout};
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

% generate a symmetric random square matrix
function M = randsym(varargin)
    M = rand(varargin{:});
    MU = triu(M);
    M = MU + MU' - diag(diag(M));
end

% force adjacency matrix to be block diagonal
function [A_block, p] = forceblockdiagonal(A)
    g = graph(sparse(A));
    bins = conncomp(g);
    [~, p] = sort(bins);
    A_block = A(p, p);
end