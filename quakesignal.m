function x = quakesignal(t, t0, r, d, F, P, A)
% x = QUAKESIGNAL(t, t0, r, d, F, P, A)
%
% Generates fake signals that looks like an earthquake signal.
%
% INPUT:
% t         time
% t0        middle onset time, "arrival time"       [default: 0]
% r         onset time scale (less = rise faster)   [default: 0.1]
% d         decay time scale (more = decay slower)  [default: 2]
% F         frequencies of the cosine waves         [default: []]
% P         phases of the cosine waves              [default: []]
% A         amplitudes of the cosine waves          [default: []]
%
% OUTPUT:
% x         signal
%                       exp(-t/d)           /                       \
%           x(t+t0) = ---------------  Sum | A_i cos( 2 pi F_i + Pi) |
%                      1 + exp(-t/r)        \                       /
%
% Last modified by spipatprathanporn@ucsd.edu, 11/26/2025

defval('t0', 0)
defval('r', 0.1)
defval('d', 2)
defval('F', [])
defval('A', [])
defval('P', [])

if size(t, 1) == 1
    t = t';
end
if size(A, 1) == 1
    A = A';
end
if size(F, 1) == 1
    F = F';
end
if size(P, 1) == 1
    P = P';
end

% envelope function
x = exp(-(t-t0) / d) ./ (1 + exp(-(t-t0) / r));
%x = exp(-(t-t0).^2 / 2 / d^2) ./ (1 + exp(-(t-t0) / r));

% cosine waves
if ~isempty(F)
    if length(A) < length(F)
        A = padarray(A, [length(F) - length(A) 0], 0, "post");
    end
    if length(P) < length(F)
        P = padarray(P, [length(F) - length(P) 0], 0, "post");
    end
    x = x .* (cos(2*pi*(t-t0) * F' + P') * A / sum(A));
end
end