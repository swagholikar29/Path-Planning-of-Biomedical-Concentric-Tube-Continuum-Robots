function [k] = inplane_bending(tubes, isCurved)
%INPLANE_BENDING uses the elastic deformation formulas found in Webster2009
%  to calculate the resulting curvature
%  INPUTS
%   tubes: array of Precurved tubes
%   isCurved: [N] optional array mapping sections that are curved
%       1 = curved, 0 = straight, -1 = skip
%  OUTPUTS
%   ks: new kappas

numTubes = length(tubes);
if ~exist('isCurved', 'var')
    isCurved = ones(1, numTubes);
end

k_num = 0;
k_dem = 0;
for i = 1:numTubes
    t = tubes(i);
    M = t.E * t.I;
    
    % skip current tube
    if isCurved(i) == -1, continue, end
    
    k_num = k_num + (M * (isCurved(i) * t.precurve));
    k_dem = k_dem + (M);
end

k = k_num/k_dem;
end

