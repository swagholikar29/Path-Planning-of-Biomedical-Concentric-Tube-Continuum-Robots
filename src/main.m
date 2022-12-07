%% Script to create series of precurved concentric tubes and evaluate their kinematics and mechanics
clear, clc, close all
addpath('path-planning/kinematics');
addpath('utils');
%% CREATE PRECURVED TUBES

% parameters for generating tubes
ODs = [4 3 2] .* 1e-3;          % (m) outer diameters 
IDs = ODs - 1.2e-3;             % (m) inner diameters
E = 1.515e9;                    % Young's Modulus

radii = [67.36796562	40.61657128];
% 
precurvature = 1000./radii;
precurvature = [15 30 45];
Ls = [75e-3 140e-3 140e-3];                    % (m) length of straight section
Lc = [30e-3 30e-3 30e-3];                     % (m) length of curved section

robot = ConcentricTubeRobot(ODs, IDs, precurvature, Ls, Lc, E);

q = [10e-3   deg2rad(0);        % outermost tube
     25e-3  deg2rad(180);       % middle tube
     40e-3      deg2rad(0)];    % innermost tube
% q = [20e-3 deg2rad(-90)];

robot.fwkine(q, false);
robot.plotTubes();
plotAllTubes(robot.tubes, q);

% r_ang = robot.arcs(:, 2, 1);

% hold on;
figure;
hold on;
h = robot.plotEnergyContour([deg2rad(88) 0]);
scatter(0.958, 1.374, 100, '*',  'k');
scatter(deg2rad(0), deg2rad(90),  100, 'x', 'k');
box('on');
Set the remaining axes properties
% 
%% Single tube
robot = ConcentricTubeRobot(3e-3, 2e-3, [30], 100e-3, 50e-3, 1.515e9);
robot.fwkine([150e-3 0], false);
robot.plotTubes();

%% Animate tube configurations
set(gcf, 'WindowState', 'maximized');
axis equal;

myVideo = VideoWriter('manipulator animated', 'MPEG-4'); %open video file
myVideo.FrameRate = 30;  %can adjust this, 5 - 10 works well for me
myVideo.Quality = 100;
open(myVideo)


n = 50;

rho2 = 25e-3;
rho2step = linspace(0, rho2, n);
theta2 = deg2rad(-90);
theta2step = linspace(0, theta2, n);

rho3 = 50e-3;
rho3step = linspace(rho2, rho3, n);
theta3 = deg2rad(90);
theta3step = linspace(0, theta3, n);

z = zeros(1, n);
one = ones(1,n);

q21 = [rho2step one.*rho2  one.*rho2   flip(rho2step)];
q22 = [z        theta2step one.*theta2 one.*theta2];
q31 = [rho2step rho3step   one.*rho3   linspace(rho3, 0, n)];
q32 = [z        z          theta3step  one.*theta3];

for i=1:length(q21)
    q = [0.0 0.0;
         q21(i) q22(i);
         q31(i) q32(i)];
    
    robot.fwkine(q, false);
    robot.animateTubes();
    
    xlim([-.0028 .0507]);
    ylim([-.03 .005]);
    zlim([0 0.0892]);
    ylabel('Y (mm)');
    xlabel('X (mm)');
    zlabel('Z (mm)');
    
    pause(0.001) %Pause and grab frame
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
end

% 
close(myVideo)
