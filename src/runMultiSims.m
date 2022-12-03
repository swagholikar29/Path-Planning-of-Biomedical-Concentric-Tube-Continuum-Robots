%% Script to automatically run a series of simulations with different parameters
%  There are several parameters that can be changed for the simulation. 
%   Simply choose the ones that should be changed and iterate over them.
%   All workspace files are saved independently, a csv with all the results
%   is added to the folder
%

%
%  Authors: Jesse F. d'Almeida <jfdalmeida@wpi.edu>
%  
% Last Version: 5/29/2020
close all, clear, clc
addpath('kinematics', 'utils', 'figure-generation', 'path-planning', ...
        'utils/stlTools/', 'utils/visibility/', 'utils/ray-casting/', '../anatomical-models');

%% Constants for the simulations
% Anatomical model definition
modelID = 'larynx1';  % ID of the anatomical model (see the `anatomical-models' folder)

% Simulation parameters
nPoints = 100; % number of configurations sampled by RRT
dq = 0.06;    % step size for rrt
laserOffsetAngle = 90; %remove
useWrist = true; %remove

% Endoscope geometry definition
%  The variable naming used in this section is consistent with (Chiluisa et al. ISMR 2020)
n = 10; % number of cutouts
viewang = deg2rad(85);
u = [1.8646e-3 * ones(1,n-1), 1.86e-3]; % notch spacing [m]
h = 0.1191e-3 * ones(1,n);         % notch height  [m]

%if useWrist
%    % wrist configuration
%    w = 0.935e-3 * ones(1,n);          % notch width   [m]
%    OD = 1.10e-3;                     % endoscope outer diameter [m]
%    ID = 0.90e-3;                     % endoscope inner diameter [m]
else
    fprintf("Wrist is turned off!!\n")
    % laser configuration
    w = 0.50e-3 * ones(1,n);          % notch width   [m]
    OD = 0.60e-3;                     % endoscope outer diameter [m]
    ID = 0.20e-3;                     % endoscope inner diameter [m]
end

%if laserOffsetAngle
%    fprintf("Laser Angle offset of %d deg\n", laserOffsetAngle)
%end

%% Parameters to vary

params = [100];

%% Set up csv to collect data
filename = 'Simulation-Larynx1-numberNotches.csv';

otherinfo = [];
if ~useWrist
    otherinfo = [otherinfo '-nowrist-'];
end

%% Loop through different sims
len = size(params, 2);
start = tic;
simIDs = {};
parfor pp = 1 : len
    nPoints = params(pp);
    
    simulationID = [modelID otherinfo '-nNotches-' num2str(n) '-' num2str(nPoints) 'pts'];    % specific workspace id
    simIDs{pp} = simulationID;
    
    fprintf('** Beginning simulation %d of %d ID: %s ** \n', pp, len, simulationID)
    calcReachableSpace(u, h, w, ID, OD, modelID, nPoints, dq, useWrist, simulationID); % runs rrt
    calcVisibleArea(simulationID, 'mcrc', laserOffsetAngle);                          % runs raycasting
    % record important data
    getSimData(simulationID, filename);
    %animateResults(simulationID);
end

%% display figures
simIDs = convertCharsToStrings(simIDs);
for i=1:length(simIDs)
    simulationID = simIDs{i};
    makeVisibilityFig(simulationID); 
    %savefig(['figures/' simulationID '.fig']);
end

fintime = toc(start);
fprintf('\n***All Simulations are complete in %.2f min. Data collected in file %s***\n', fintime/60, filename)
