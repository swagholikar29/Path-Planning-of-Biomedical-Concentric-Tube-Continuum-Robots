%% Script to create series of precurved concentric tubes and evaluate their kinematics and mechanics
close all, clear, clc
addpath('kinematics');

%% CREATE PRECURVED TUBES
nTubes = 2;
ODs = [5e-3 3e-3 2e-3];
IDs = ODs - 1e-3;
precurves = [15 30 90];
Ls = 30e-3;
Lc = 30e-3;

for i = 1:nTubes
    tubes(i) = Precurved(ODs(i), IDs(i), precurves(i), Ls, Lc);
    tubes(i).fwkine(deformation(tubes(i).precurve, Lc, 10e-3));
end

% plot the tubes individually
plotTubes(tubes);

% deformed together 
acts = [10e-3 0; 
        30e-3 0];
q = actuator2arcparams(tubes, acts)
for i = 1:nTubes
    tubes(i).fwkine(q(:,:,i));
end

%% PLOT
plotTubes(tubes);


