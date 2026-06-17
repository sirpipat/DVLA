function [g, A] = fakesnowflakes(n, s, kmin)
% [g, A] = fakesnowflakes(n, s, kmin)
%
% Randomly generates a graph made of smaller separated subgraphs called
% snowflakes.
%
% INPUT:
% n             number of nodes in the graph
% s             number of subgraph
% kmin          minimum number of nodes in a subgraph
%
% OUTPUT:
% g             graph
% A             connectivity matrix
%
% Last modified by spipatprathanporn@ucsd.edu, 06/09/2026

defval('kmin', 2)

% randomly assign the subgroup index to nodes. (s*kmins nodes have been 
% taken out, so that every subgroup is guaranteed to have at least kmin 
% nodes.)
ni = randi(s, n-s*kmin, 1);

% group nodes with their assigned indices into subgroups
ns = histcounts(ni, 'BinEdges', (0:s)+0.5) + kmin;

A = zeros(n);
for ii = 1:length(ns)
    jj_begin = sum(ns(1:(ii-1))) + 1;
    jj_end = sum(ns(1:ii));
    A(jj_begin:jj_end, jj_begin:jj_end) = ones(ns(ii)) - eye(ns(ii));
end
g = graph(A);
end