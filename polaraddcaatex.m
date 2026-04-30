function sc = polaraddcaatex(ax)
% sc = POLARADDCATEX(ax)
%
% Adds the locations of the CAATEX hydrophone moorings to a polar map.
%
% INPUT:
% ax        target axes: where you want to add 
%
% OUTPUT:
% sc        object handles to the scatter plots
%
% SEE ALSO:
% POLARBASEMAP, PLOTPOLARMAP
%
% Last modified by spipatprathanporn@ucsd.edu, 04/30/2026

defval('ax', gca)

% mooringinfo
kstnm = {'SIO1', 'SIO2', 'SIO3', 'NERSC1', 'NERSC2', 'NERSC3'}';
stla = [80.8697, 77.7301, 73.0180, 84.0033, 82.5125, 83.4420]';
stlo = [-148.3434, -149.1429, -149.5776, 28.3916, 23.9628, 25.6672]';
stco = [1 0.6 0.2; 0 1 1; 1 0.5 1; 1 0 0; 0 1 0; 1 1 0];

% create object placeholder matrix
sc = gobjects(6,1);

% scatter plot
for ii = 1:6
    sc(ii) = plotpolarmap(@scatter, ax, stlo(ii), stla(ii), 100, ...
        stco(ii,:), 'filled', 'v', 'MarkerEdgeColor', 'k', ...
        'DisplayName', kstnm{ii});
end
end