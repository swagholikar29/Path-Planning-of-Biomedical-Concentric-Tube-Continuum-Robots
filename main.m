%% Script to create series of precurved concentric tubes and evaluate their kinematics and mechanics

close all, clear, clc
addpath('kinematics');

outer = Precurved(5e-3, 15, 30e-3, 30e-3);
q = [0 10e-3];
outer.fwkine(q);

inner = Precurved(2e-3, 50, 30e-3, 30e-3);
q = [0 30e-3];
inner.fwkine(q);

tubes = [outer inner];

plotTubes(tubes);


