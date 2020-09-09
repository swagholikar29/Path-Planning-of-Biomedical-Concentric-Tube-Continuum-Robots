%% Script to create series of precurved concentric tubes and evaluate their kinematics and mechanics
close all, clear, clc
addpath('kinematics');

%% CREATE PRECURVED TUBES

outer = Precurved(5e-3, 4e-3, 15, 30e-3, 30e-3);
q1 = [0 10e-3]; 

inner = Precurved(2e-3, 1e-3, 50, 30e-3, 30e-3);
q2 = [0 10e-3];
tubes = [outer inner];

k = inplane_bending(tubes)

q1 = [k q1 0];
q2 = [k q2 0];
outer.fwkine(q1);
inner.fwkine(q2);

%% PLOT
plotTubes(tubes);


