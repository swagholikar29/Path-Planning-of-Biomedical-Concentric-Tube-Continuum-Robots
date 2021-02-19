%% Script to create series of precurved concentric tubes and evaluate their kinematics and mechanics
close all, clear, clc
addpath('kinematics');

%% CREATE PRECURVED TUBES

% parameters for generating tubes
ODs = [5.4e-3 3.8e-3];         % (m) outer diameters 
IDs = ODs - 0.6e-3;             % (m) inner diameters
E = 1900e6;

precurves = [30 15];
Ls = [100-3 250-3];                    % (m) length of straight section
Lc = [50e-3 50e-3];                     % (m) length of curved section

robot = ConcentricTubeRobot(ODs, IDs, precurves, Ls, Lc, E);

q = [20e-3   deg2rad(90);        % outermost tube
     20e-3  deg2rad(0)];     % innermost tube
 
robot.fwkine(q, false);
robot.plotTubes();
arcs = robot.arcs;
disp('Curvature (1/m)');
disp(arcs(:,1,end));

disp('Rotation (deg)');
disp(rad2deg(arcs(:,2,end)));


