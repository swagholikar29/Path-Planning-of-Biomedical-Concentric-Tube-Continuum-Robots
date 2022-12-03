%% Script to test the kinematics of the the Endoscope
% clc, clear all
% addpath('kinematics')
% addpath('path-planning')
% addpath('utils')

%% Load in Simulation and STL
simulationID = 'larynx1-TEST--Laser_ang-90-dq-0.06-100pts.mat';
load(simulationID);

% load in cropped version of stl mesh
path = fullfile('..', 'anatomical-models', modelID);
pathStl = fullfile('..', 'anatomical-models', modelID, 'tissue_cropped.stl');
[vertices, faces, ~, ~] = stlRead(pathStl);
earModel.vertices = vertices .* 1e-3;
earModel.faces = faces;

%% Fwkin for specific index

idx = 470;
lim = 75;

laserAng = deg2rad(90);           % convert laser into degrees
laserVec = [sin(laserAng) 0 cos(laserAng) 0]';  % unit vector of laser direction (0 at end so that it can be transformed)

total_visibility = [];
for ii = lim 
    p = pList(:,ii);
    q = qList(:,ii);
    robot.fwkine(qList(:,ii), earModel.baseTransform);
    
    T_tip = TList(:,:,ii);
    approach4 = T_tip * laserVec;                   % approach vector with extra element from transform
    approach = approach4(1:3);                      % split 3D vector

    % Creat Visibility Map
%     approach = robot.transformations(1:3,3,end);    % z direction of transform to tip
    [vmap, quiver] = visibilitymap(p, approach, earModel, 'mcrc');
    total_visibility(:, ii) = vmap;
    
end

total_visibility = logical(sum(total_visibility, 2));

%% Plot model
figure
hold on; grid on; axis equal;

robotPhysicalModel = robot.makePhysicalModel();

% STL plot
h1 = stlPlot(earModel.vertices, earModel.faces, '', total_visibility);

% Wrist plot
h2 = surf(robotPhysicalModel.surface.Xw, ...
    robotPhysicalModel.surface.Yw, ...
    robotPhysicalModel.surface.Zw, ...
    'FaceColor','blue');

% Endoscope plot
h3 = surf(robotPhysicalModel.surface.Xe, ...
    robotPhysicalModel.surface.Ye, ...
    robotPhysicalModel.surface.Ze, ...
    'FaceColor','red');

% tip position plot
scatter3(pList(1,lim), pList(2,lim), pList(3,lim), 'filled', 'red');

vp = quiver.vp * 1e-3;
rays = quiver.dispR * 1e-3;
q1 = quiver3(vp(1,:), vp(2,:), vp(3,:), rays(1,:), rays(2,:), rays(3,:));
q1.Color = 'c';
    
xlabel('X [m]');
ylabel('Y [m]');
zlabel('Z [m]');
title([num2str(idx) 'th Pose'])
set(gca,'FontSize',12);

