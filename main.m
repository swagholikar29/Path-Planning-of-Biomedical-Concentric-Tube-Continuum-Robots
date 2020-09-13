%% Script to create series of precurved concentric tubes and evaluate their kinematics and mechanics
close all, clear, clc
addpath('kinematics');

%% CREATE PRECURVED TUBES
nTubes = 3;
ODs = [5e-3 4e-3 3e-3 2e-3];
IDs = ODs - 1e-3;
precurves = [15 50 70 90];
Ls = 100e-3;
Lc = 30e-3;

% make the precurved tubes
for i = 1:nTubes
    tubes(i) = Precurved(ODs(i), IDs(i), precurves(i), Ls, Lc);
end

% Joint Parameters of the tube
%  p: trans  alpha: rotation
q = [10e-3   deg2rad(0); 
     30e-3   deg2rad(60);
     50e-3   deg2rad(90)
     40e-3   deg2rad(0)];

arcs = joint2arcparams(tubes, q);
for i = 1:nTubes
    tubes(i).fwkine(arcs(:,:,i));
end

%% PLOT
plotTubes(tubes);
plotAllTubes(tubes, q);


