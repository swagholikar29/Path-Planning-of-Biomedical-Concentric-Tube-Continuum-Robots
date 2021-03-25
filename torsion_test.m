
%% CREATE PRECURVED TUBES
clear, clc, close all;
% parameters for generating tubes
ODs = [2.39e-3 1.6e-3];         % (m) outer diameters 
IDs = [2.01e-3 0];             % (m) inner diameters
E = 75e6;

precurves = [.0099e3 .0138e3];
Ls = [93.5e-3 218.5e-3];                    % (m) length of straight section
Lc = [92.3e-3 85e-3];                     % (m) length of curved section

robot = ConcentricTubeRobot(ODs, IDs, precurves, Ls, Lc, E);
robot.v = .35;

q = [10e-3   deg2rad(88);        % outermost tube
     10e-3  deg2rad(0)];     % innermost tube

robot.fwkine(q, false);
% robot.plotTubes();

figure
degs = 90:45:270;

for i=1:length(degs)
    
    subplot(round(length(degs)/2), 2, i);
    
    robot.plotEnergyContour([deg2rad(0) deg2rad(degs(i))]);
end
