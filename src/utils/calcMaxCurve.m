function kappa = calcMaxCurve(D,maxStrain)
%CALCMAXCURVE Uses yield strength and tube diameter
%to calculate the precurvature manufactured into tube
%INPUTS:
%   D = outer diameter of tube
%   maxStrain = strain at max yield strength for material

kappa = 2*maxStrain / (D * (1 + maxStrain));
end

