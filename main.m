%% Script to create series of precurved concentric tubes and evaluate their kinematics and mechanics
close all, clear, clc
addpath('kinematics');

%% CREATE PRECURVED TUBES

% parameters for generating tubes
nTubes = 2;                     % number of tubes
ODs = [3.25e-3 2e-3];         % (m) outer diameters 
IDs = [2.45e-3 1.2e-3];               % (m) inner diameters
precurves = [30 50];         % (1/m) curvatures
Ls = 50-3;                    % (m) length of straight section
Lc = 50e-3;                     % (m) length of curved section

% make the precurved tubes
for i = 1:nTubes
    tubes(i) = Precurved(ODs(i), IDs(i), precurves(i), Ls, Lc);
end

%% Calculate tube deformations and forward kinematics

% -------------change these parameters-------------
% Joint Parameters of the tube
% translation  rotation
q = [0e-3   deg2rad(0);        % outermost tube
     20e-3   deg2rad(30)];     % innermost tube

% calculate forward kinematics
arcs = joint2arcparams(tubes, q)   % calcs arc parameters of the deformation
for i = 1:nTubes
    tubes(i).fwkine(arcs(:,:,i));
end

%% PLOT
plotTubes(tubes);           % regular plot of deformed tubes
plotAllTubes(tubes, q);     % uncomment for plot of each tube undeformed


