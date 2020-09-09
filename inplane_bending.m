function [k] = inplane_bending(tubes)
%INPLANE_BENDING uses the elastic deformation formulas found in Webster2009
%  to calculate the resulting curvature
%  INPUTS
%   tubes: array of Precurved tubes
%   qs: matrix of q for each tube [base_rotation dz]
%  OUTPUTS
%   ks: new kappas

numTubes = length(tubes);

ks = zeros(numTubes, 1);

k = 0;

for i = 1:numTubes
    t = tubes(i);
    newK = (t.E*t.I*t.precurve)/(t.E*t.I);
    k = k + newK;
end
end

