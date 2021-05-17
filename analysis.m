%% Calc Predicted Value
clear, clc, close all
addpath('kinematics');

%% Create manipulator obj

% parameters for generating tubes
ODs = [5.6 3.9].* 1e-3;         % (m) outer diameters 
IDs = ODs - 1.2e-3;             % (m) inner diameters
E = 1.515e9;

radii = [34.22853397	62.01856163];
precurves = 1000./radii;
Ls = [75e-3 140e-3];                    % (m) length of straight section
Lc = [50e-3 50e-3];                     % (m) length of curved section

robot = ConcentricTubeRobot(ODs, IDs, precurves, Ls, Lc, E);

%% FW kin
% get initial tip position
q = [0e-3   deg2rad(88);        % outermost tube
     0e-3  deg2rad(0)];       % innermost tube
robot.fwkine(q, false);
% robot.plotTubes();
disp(rad2deg(robot.arcs(3,2,1)))
% hold on;

robot.fwkine(q, true);
% robot.plotTubes();
disp(rad2deg(robot.arcs(3,2,1)))
%%
% tip0 = robot.tubes(end).pose(:,end);
% 
% % move robot
% q = [0e-3   deg2rad(0);        % outermost tube
%      20e-3  deg2rad(0)];       % innermost tube
% 
% robot.fwkine(q, false);
% % robot.plotTubes();
% 
% % get end tip position
% tip = robot.tubes(end).pose(:,end);
% R2world = [0  0 -1; 
%            -1  0 0; 
%            0 -1 0];         % rotation M from robot to world
% 
% tip_diff = tip - tip0;
% tip_world = R2world*tip_diff*1e3;
% disp(tip_world);


