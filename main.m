%% Script to create series of precurved concentric tubes and evaluate their kinematics and mechanics
close all, clear, clc
addpath('kinematics');

%% CREATE PRECURVED TUBES
ODs = [5e-3 2e-3];
IDs = ODs - 1e-3;
precurves = [15 50];
Ls = 30e-3;
Lc = 30e-3;

outer = Precurved(ODs(1), IDs(1), precurves(1), Ls, Lc);
inner = Precurved(ODs(2), IDs(2), precurves(1), Ls, Lc);
tubes = [outer inner];

k = inplane_bending(tubes)

outer.fwkine(outer.deformation(k, Lc, 10e-3));
inner.fwkine(inner.deformation(k, Lc, 10e-3));

%% PLOT
plotTubes(tubes);


