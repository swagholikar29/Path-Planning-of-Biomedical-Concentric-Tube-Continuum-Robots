%% Script to create series of precurved concentric tubes and evaluate their kinematics and mechanics
close all, clear, clc
addpath('kinematics');

%% CREATE PRECURVED TUBES

% parameters for generating tubes
ODs = [5.6e-3 3.9e-3];         % (m) outer diameters 
IDs = ODs - 0.6e-3;             % (m) inner diameters
E = 1.515e9;

radii = [33.3 66.7];

precurves = 1000./radii;
Ls = [100-3 250-3];                    % (m) length of straight section
Lc = [50e-3 50e-3];                     % (m) length of curved section

robot = ConcentricTubeRobot(ODs, IDs, precurves, Ls, Lc, E);

q = [20e-3   deg2rad(88);        % outermost tube
     20e-3  deg2rad(0)];     % innermost tube

robot.fwkine(q, false);
robot.plotTubes();

angles = zeros(180, 1);
pause(.1);

% for i = 1:180
%     
%     q = [20e-3   deg2rad(i);        % outermost tube
%          20e-3  deg2rad(0)];     % innermost tube
%     robot.fwkine(q, false);
%     robot.animateTubes();
%     
%     arcs = robot.arcs;
%     angles(i) = rad2deg(arcs(3,2,end));
%     
%     disp(['Input ' num2str(i) ' Output ' num2str(angles(i))])
% %     pause(.1);
% end
plot(1:180, angles)
xlabel('Input Angle (deg)')
ylabel('Output Angle (deg)')
title('Angle of tubes with torision')

arcs = robot.arcs;
disp('Curvature (1/m)');
disp(arcs(:,1,end));

disp('Rotation (deg)');
disp(rad2deg(arcs(:,2,end)));


