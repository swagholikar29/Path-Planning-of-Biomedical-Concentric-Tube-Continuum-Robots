%% Script to create series of precurved concentric tubes and evaluate their kinematics and mechanics
close all, clear, clc
addpath('kinematics');

%% CREATE PRECURVED TUBES

% parameters for generating tubes
nTubes = 2;                     % number of tubes
ODs = [6e-3 4e-3];         % (m) outer diameters 
IDs = ODs - 1.6e-3;
% IDs = [2.4e-3 0.4e-3];               % (m) inner diameters

radii = [81.27];               
k = 1 ./ (radii * 1e-3);         % (1/m) curvatures

precurves = [0 k];

Ls = 100-3;                    % (m) length of straight section
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
     0e-3  deg2rad(0)];     % innermost tube

% calculate forward kinematics
arcs = joint2arcparams(tubes, q);   % calcs arc parameters of the deformation
fprintf("Resultant Radius of Curvature: %.2fmm \n", 1/(arcs(3, 1, 1) *1e-3));
for i = 1:nTubes
    tubes(i).fwkine(arcs(:,:,i));
end

%% PLOT
% plotTubes(tubes);           % regular plot of deformed tubes
% plotAllTubes(tubes, q);     % uncomment for plot of each tube undeformed


