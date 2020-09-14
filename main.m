%% Script to create series of precurved concentric tubes and evaluate their kinematics and mechanics
close all, clear, clc
addpath('kinematics');

%% CREATE PRECURVED TUBES

% parameters for generating tubes
nTubes = 3;                     % number of tubes
ODs = [5e-3 4e-3 3e-3];         % (m) outer diameters 
IDs = ODs - 1e-3;               % (m) inner diameters
precurves = [15 50 70];         % (1/m) curvatures
Ls = 100e-3;                    % (m) length of straight section
Lc = 30e-3;                     % (m) length of curved section

% make the precurved tubes
for i = 1:nTubes
    tubes(i) = Precurved(ODs(i), IDs(i), precurves(i), Ls, Lc);
end

%% Calculate tube deformations and forward kinematics

% -------------change these parameters-------------
% Joint Parameters of the tube
% translation  rotation
q = [0e-3   deg2rad(0);        % outermost tube
     10e-3   deg2rad(30);
     20e-3   deg2rad(90)];     % innermost tube

% calculate forward kinematics
arcs = joint2arcparams(tubes, q);   % calcs arc parameters of the deformation
for i = 1:nTubes
    tubes(i).fwkine(arcs(:,:,i));
end

%% PLOT
plotTubes(tubes);           % regular plot of deformed tubes
plotAllTubes(tubes, q);     % uncomment for plot of each tube undeformed


