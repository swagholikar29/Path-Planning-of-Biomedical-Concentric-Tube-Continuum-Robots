%% Script to characterize workspace of the Endoscope or EndoWrist
clc, clear, close all
addpath('../kinematics')
addpath('../path-planning')
addpath('../utils')

% "Endoscope", "EndoWrist" or "Wrist"
robotType = "EndoWrist";

nPoints = 75e3;
delta_q = [0.009 0.02 0.00275 0.0175 0.015 0.00275];

xpos = [];
ypos = [];
zpos = [];

% Parameters for EndoWrist
n = 5; % number of cutouts
u = [0.0025 * ones(1,n-1), 2e-3]; % notch spacing [m]
h = 4.7172e-04 * ones(1,n);       % notch height  [m]
w = 1.40e-3 * ones(1,n);          % notch width   [m]
OD = 1.60e-3;                     % endoscope outer diameter [m]
ID = 1.40e-3;                     % endoscope inner diameter [m]
cutouts.w = w;
cutouts.u = u; % converting to meters
cutouts.h = h; % converting to meters
cutouts.alpha = zeros(1,length(u));

maxBend = 0.025; % [1/m]
maxKappa= 1/maxBend;
maxTheta= 4*pi; % [rad]
maxDz   = 20e-3; %[m]
maxDisp = sum(cutouts.h); % [m]
maxRot  = 4*pi;  % [rad]
maxAdv  = 20e-3; % [m]

qBounds = [0        0        0     0       0       0;
           maxKappa maxTheta maxDz maxDisp maxRot  maxAdv];

if robotType == "Endoscope"
    robot = Endoscope();
    delta_q = delta_q(1:3);
    qBounds = qBounds(:,1:3);
    qList = new_rrt(nPoints, delta_q, qBounds);
    disp('Calculating Poses')
    for n = 1:nPoints
        robot.fwkine(qList(:,n));
        xpos = [xpos robot.pose(1, end)];
        ypos = [ypos robot.pose(2, end)];
        zpos = [zpos robot.pose(3, end)];
    end
elseif robotType == "EndoWrist" % robotType is EndoWrist
    robot = EndoWrist(ID, OD, n, cutouts);
    qList = new_rrt(nPoints, delta_q, qBounds);
    disp('Calculating Poses')
    for n = 1:nPoints
        robot.fwkine(qList(:,n));
        xpos = [xpos robot.pose(1)];
        ypos = [ypos robot.pose(2)];
        zpos = [zpos robot.pose(3)];
    end
elseif robotType == "Wrist"
    robot = Wrist(ID, OD, n, cutouts);
    delta_q = delta_q(4:6);
    qBounds = qBounds(:,4:6);
    qList = new_rrt(nPoints, delta_q, qBounds);
%    [qListNormalized,qList,pList,aList] = rrt_ndim(robot, qBounds, [], delta_q, nPoints);
    disp('Calculating Poses')
    for n = 1:nPoints
        robot.fwkine(qList(:,n), eye(4));
        xpos = [xpos robot.pose(1, end)];
        ypos = [ypos robot.pose(2, end)];
        zpos = [zpos robot.pose(3, end)];
    end
end
 
% disp('Making Scatterplot')
% figure
% scatter3(xpos, ypos, zpos, '.', 'MarkerFaceColor', '#9400D3', 'MarkerEdgeColor', '#9400D3')
% title([robotType, " Workspace - ", num2str(nPoints), " Points"])
% xlabel('X [m]')
% ylabel('Y [m]')
% zlabel('Z [m]')
% axis equal

figure
plot(xpos, zpos, '.', 'MarkerEdgeColor', '#9400D3')
xlabel('X [m]')
zlabel('Z [m]')
title([robotType, " Cutout View - ", num2str(nPoints), " Points"])
axis equal