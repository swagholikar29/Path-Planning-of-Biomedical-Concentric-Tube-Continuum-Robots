%% This script simulates the exploration with a endoscope-wrist robot.
%  The algorithm works as follows: first, we run a sampling-based path planning
%  algorithm (RRT) to generate a set of reachable points within space
% (based off runme)
%
%  Authors: A. Chiluisa <ajchiluisa@wpi.edu>
%           L. Fichera  <lfichera@wpi.edu>
%           I. Chan <iachan@wpi.edu>
%  
% Last Version: 5/29/2020
close all, clear, clc
addpath('kinematics', 'utils', 'figure-generation', 'path-planning', ...
        'utils/stlTools/', 'utils/visibility/', 'utils/ray-casting/', ...
        'utils/wrist_configs/', '../anatomical-models', 'simAnalyzer/');
    
%% Simulation parameters
nPoints = 10; % number of configurations sampled by RRT
dq = 0.06;
useWrist = true; %To be removed
laserOffsetAngle = 0; %To be removed

%% Anatomical model definition
modelID = 'larynx1'; % ID of the anatomical model (see the `anatomical-models' folder)
%%%To be removed%%%
otherinfo = [];
if ~useWrist %To be removed
    otherinfo = [otherinfo '-nowrist']; 
end
if laserOffsetAngle
    otherinfo = [otherinfo '-Laser_ang-' num2str(laserOffsetAngle) '-'];
end
%%%


simulationID = [modelID otherinfo '-dq-' num2str(dq) '-' num2str(nPoints) 'pts'];
simulationID = [modelID otherinfo '-dq-' num2str(dq) '-' num2str(nPoints) 'pts'];



%% Endoscope geometry definition
%  The variable naming used in this section is consistent with (Chiluisa et al. ISMR 2020)
n = 10; % number of cutouts
viewang = deg2rad(85);
R = 8;
%L = calc_L(viewang, R, 2.5, 0, 0);
L = 15;
[singleH, singleU] = calc_config(L, R, n, 1.1/2, 0.9/2, 0.935);

%u = 0.000367766480556499;
%h = 0.000187789074999056;
u = singleU*1e-3 * ones(1,n); % notch spacing [m]
h = singleH*1e-3 * ones(1,n);         % notch height  [m]

if useWrist
    % wrist configuration
    w = 0.935e-3 * ones(1,n);          % notch width   [m]
    OD = 1.1e-3;                     % endoscope outer diameter [m]
    ID = 0.9e-3;                     % endoscope inner diameter [m]
else
    fprintf("Wrist is turned off!\n")
    % laser configuration
    w = 0.40e-3 * ones(1,n);          % notch width   [m]
    OD = 0.60e-3;                     % endoscope outer diameter [m]
    ID = 0.40e-3;                     % endoscope inner diameter [m]
end

if laserOffsetAngle
    fprintf("Laser Angle offset of %d deg\n", laserOffsetAngle)
end

%% Estimate the reachable workspace with RRT`
calcReachableSpace(u, h, w, ID, OD, modelID, nPoints, dq, useWrist, simulationID);

%% Remove Points
RemovePoints(simulationID);

%% Create a video of this simulation and Histogram
makeVisibilityFig(simulationID);

filename = 'testSim.csv';
getSimData(simulationID, filename);

% Save figure to folder
savefig(['figures/' simulationID '.fig']);

animateResults(simulationID); 
