%% Script to create series of precurved concentric tubes and evaluate their kinematics and mechanics
close all, clear, clc
addpath('kinematics');

%% CREATE PRECURVED TUBES

% parameters for generating tubes
ODs = [3e-3 2e-3];         % (m) outer diameters 
IDs = ODs - 0.5e-3;             % (m) inner diameters
E = 1900e6;

precurves = [15 30];
Ls = [100e-3 100e-3];                    % (m) length of straight section
Lc = [50e-3 50e-3];                     % (m) length of curved section

robot = ConcentricTubeRobot(ODs, IDs, precurves, Ls, Lc, E);

q = [0e-3   deg2rad(0);        % outermost tube
     20e-3  deg2rad(180)];     % innermost tube
 
robot.fwkine(q);
robot.plotTubes();


