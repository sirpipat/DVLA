function C = splitnan(M)
% C = SPLITNAN(M)
%
% Split the matrix M at all-NaN rows. See
% https://www.mathworks.com/matlabcentral/answers/
% 335333-split-array-into-separate-arrays-at-row-with-nans
% for the method.
%
% INPUT:
% M             Matrix
%
% OUTPUT:
% C             cell array containing splitted matrix
%
% EXAMPLE:
% M = [1 2 3; 4 5 6; NaN NaN NaN; 7 8 9; 10 NaN 12; 13 14 15];
% C = splitnan(M);
%
% % Expected result
% % C = {[1 2 3; 4 5 6], [7 8 9; 10 NaN 12; 13 14 15})
%
% SEE ALSO:
% ACCUMARRAY, SPLIT
%
% Last modified by spipatprathanporn@ucsd.edu, 03/19/2026

idx = all(isnan(M), 2);
idy = 1 + cumsum(idx);
idz = 1:size(M, 1);
C = accumarray(idy(~idx), idz(~idx), [], @(r) {M(r,:)});
end