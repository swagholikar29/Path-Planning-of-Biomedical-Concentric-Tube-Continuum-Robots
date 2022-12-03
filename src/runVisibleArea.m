clear all; clc; close all;

%% Load the mat file that you want to run

load('larynx1-nowrist-dq-0.06-150pts.mat');

%%
% First, load the mat file

%% offste angle [0 45 70 90] deg
useWrist = false;
laserOffsetAngle = 25;

%% Anatomical model definition

otherinfo = [];

if ~useWrist
    otherinfo = [otherinfo '-nowrist-'];
end

if laserOffsetAngle
    otherinfo = [otherinfo 'Laser_ang-' num2str(laserOffsetAngle) '-'];
end

simulationID = [modelID otherinfo 'dq-' num2str(dq) '-' num2str(nPoints) 'pts'];

if laserOffsetAngle
    fprintf("Laser Angle offset of %d deg\n", laserOffsetAngle)
end

%% params to keep track of
results.modelID = modelID;
results.nPoints = nPoints;
results.dq = dq;
results.simID = simulationID;
%visibleMap = [];
save([simulationID '.mat'], '-v7.3');

%% Run Ray casting
tstart = tic;
calcVisibleArea(simulationID, 'mcrc', laserOffsetAngle);
telapsed = toc(tstart);
s = seconds(telapsed);
s.Format = 'hh:mm:ss'
%% Create a video of this simulation and Histogram
makeVisibilityFig(simulationID);
filename = 'testSim.csv';
getSimData(simulationID, filename);
% Save figure to folder
savefig(['figures/' simulationID '.fig']);
animateResults(simulationID);