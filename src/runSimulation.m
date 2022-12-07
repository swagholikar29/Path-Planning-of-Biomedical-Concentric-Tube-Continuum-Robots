%% This script simulates the exploration with a endoscope-wrist robot.
%  The algorithm works as follows: first, we run a sampling-based path planning
%  algorithm (RRT) to generate a set of reachable points within space
% (based off runme)

close all, clear, clc
addpath('kinematics', 'utils', 'figure-generation', 'path-planning', ...
        'utils/stlTools/', 'utils/visibility/', 'utils/ray-casting/', ...
        'utils/wrist_configs/', '../anatomical-models', 'simAnalyzer/');
    
%% Simulation parameters
nPoints = 10; % number of configurations sampled by RRT
dq = 0.06;

%% Anatomical model definition
modelID = 'larynx1'; % ID of the anatomical model (see the `anatomical-models' folder)
simulationID = [modelID '-dq-' num2str(dq) '-' num2str(nPoints) 'pts'];%%%

%% Define the robot



%% Estimate the reachable workspace with RRT`
calcReachableSpace(OD, ID, k, Ls, Lc, E, modelID, nPoints, dq, simulationID)

%% Remove Points
RemovePoints(simulationID);

%% Create a video of this simulation and Histogram
filename = 'testSim.csv';
getSimData(simulationID, filename);

% Save figure to folder
savefig(['figures/' simulationID '.fig']);

animateResults(simulationID); 
