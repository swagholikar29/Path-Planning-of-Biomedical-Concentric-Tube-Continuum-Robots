function kappa = calcMaxCurve(D,maxStrain)
%UNTITLED Summary of this function goes here
%INPUTS:
%   D = outer diameter of tube
%   maxStrain = strain at max yeild strength for material

kappa = 2*maxStrain / (D * (1 + maxStrain));
end

