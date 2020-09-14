function [k, phi] = inplane_bending(tubes, theta, isCurved)
%INPLANE_BENDING uses the elastic deformation formulas found in Webster2009
%  to calculate the resulting curvature
%  INPUTS
%   tubes    = [N] array of Precurved tubes
%   thetas   = [N] array of base rotations
%   isCurved = [N] optional array mapping sections that are curved
%       1 = curved, 0 = straight, -1 = skip
%  OUTPUTS
%   k   =  curvature of link
%   phi =  rotation of link

numTubes = length(tubes);

if ~exist('theta', 'var')
    theta = zeros(1, numTubes);
end

if ~exist('isCurved', 'var')
    isCurved = ones(1, numTubes);
end

kx_num = 0;
kx_dem = 0;

ky_num = 0;
ky_dem = 0;
for i = 1:numTubes
    t = tubes(i);
    M = t.E * t.I;
    
    % skip current tube
    if isCurved(i) == -1, continue, end
    
    kx_num = kx_num + (M * (isCurved(i) * t.precurve) * cos(theta(i)));
    kx_dem = kx_dem + (M);
    
    ky_num = ky_num + (M * (isCurved(i) * t.precurve) * sin(theta(i)));
    ky_dem = ky_dem + (M);
end

kx = kx_num/kx_dem;
ky = ky_num/ky_dem;

k = sqrt(kx.^2 + ky.^2);
phi = atan2(ky, kx);
end
