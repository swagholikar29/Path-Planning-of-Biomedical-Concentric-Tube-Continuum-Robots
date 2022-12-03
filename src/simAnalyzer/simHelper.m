%% SIMHELPER script to assemble data, figures, and videos from already ran simulation files
% Put all simulation files into one folder and change script to choose
% which function is to be ran
% 
%
% Author:   Jesse F. d'Almeida <jfdalmeida@wpi.edu>
%         
% Last revision: 7/1/2020
% 
clear all; close all; clc
simStruct = dir('simAnalyzer/simulations/*.mat'); % must run while in src folder
simIDs = convertCharsToStrings({simStruct.name});
simIDs = erase(simIDs, '.mat');
% filename = 'Trial11-Larynx1-NewOre.csv';

for i=1:length(simIDs)
    simulationID = char(simIDs(i));
    fprintf('Running on Simulation: %s\n', simulationID)
    
    animateResults(simulationID);
%     calcVisibleArea(simulationID, 'mcrc', false);
%     makeVisibilityFig(simulationID);
%     getSimData(simulationID, filename);
%     savefig(['figures/' simulationID '.fig']); 
end

