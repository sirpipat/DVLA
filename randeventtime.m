function t_event = randeventtime(t0, t1, tau, p, n, b)
% t_event = randeventtime(t0, t1, tau, p, n, b)
%
% Randomly generates the random event time series by generating the
% mainshocks first and then aftershocks. The mainshocks follow exponential
% distribution. The aftershocks follow the Omori-Utsu's law with the
% strength scaled with the time between the lastest mainshock and the
% previous one. Note that this does not generate the entire earthquake
% sequences but rather samples events from the underlying distribution.
% 
% INPUT:
% t0        beginning time of the time series           [Default: 0]
% t1        end time of the time series                 [Default: 1460.97]
% tau       recurrence time of the mainshocks           [Default: 730.985]
% p         exponent in the Omori-Utsu's law            [Default: 1]
% n         number of sampled events                    [Default: 1]
% b         intervals to exclude (mask out). It is a monotonically
%           increasing vector in this form [b1_begin b1_end b2_begin b2_end
%           ... bN_begin bN_end]                        [Default: []]
%
% OUTPUT:
% t_event   unsorted randomly generated event origin times of length n
% 
% Last modified by spipatprathanporn@ucsd.edu, 06/17/2026

defval('t0', 0)
defval('t1', t0 + 4*365.2425)   % 4 years
defval('tau', 2*365.2425)       % 2 years
defval('p', 1)
defval('n', 1)
defval('b', [])

% the beginning of time
t00 = t1 - 2*(t1 - t0);

% mainshock intervals
tau_main = -tau * log(1 - rand());
cum_tau_main = t00 + tau_main;
num_mains = 1;

% keep generating mainshock intervals until the cumulative exceed the
% duration
while cum_tau_main(end) < t1 - t00
    num_mains = num_mains + 1;
    tau_main(num_mains) = -tau * log(1 - rand());
    cum_tau_main(num_mains) = cum_tau_main(num_mains-1) + tau_main(num_mains);
end
% remove the last mainshock since it occurs after the end
num_mains = num_mains - 1;
tau_main = tau_main(1:num_mains);
cum_tau_main = cum_tau_main(1:num_mains);

% compute energy released
E_before = zeros(num_mains+1,1);
E_after = zeros(num_mains+1,1);
K = zeros(num_mains,1);

for ii = 1:num_mains
    E_before(ii+1) = E_after(ii) + 1/tau * tau_main(ii);
    % energy released: scaled to the stored energy
    % This random function is arbitrary. TODO: find a realistic model
    K(ii) = (0.7 + 0.3*rand()) * E_before(ii+1);
    E_after(ii+1) = E_before(ii+1) - K(ii);
end

% if all events happens before t0, move the latest event forward to t0
if cum_tau_main(end) < t0
    tau_main(1) = tau_main(1) + (t0 - cum_tau_main(end));
    cum_tau_main = cum_tau_main + (t0 - cum_tau_main(end));
end

% time samples
t = (floor(t00):ceil(t1))';

% the comb function: f(t) = K_i * delta(t - t_i), t_i = Sum(j=1,i) tau_j
f = 0*t;
for ii = 1:num_mains
    f(t == round(cum_tau_main(ii))) = K(ii);
end

% the exponential function decay
t_expo = (1:length(t))';
f_expo = t_expo .^ -p;

% convolve
f = conv(f, f_expo);
f = f(1:length(t));

% apply boxcar window to exclude specific periods
% This is the unnormalized probability distribution function (pdf)
for ii = 1:2:length(b)-1
    f(and(t > b(ii), t < b(ii+1))) = 0;
end

% compute the cumulative distribution function (cdf)
cdf = cumtrapz(t(and(t>=t0, t<=t1)), f(and(t>=t0, t<=t1)));
% normalize
cdf = cdf / cdf(end);

% add a tiny value to resolve repeating values in the cdf
epsilon = 1e-12;
for ii = 2:length(cdf)
    if cdf(ii) <= cdf(ii-1)
        cdf(ii) = cdf(ii-1) + epsilon;
    end
end
t_used = t(and(t>=t0, t<=t1));
% generate random event time
t_event = interp1(cdf, t_used, rand(n, 1), 'linear');
end