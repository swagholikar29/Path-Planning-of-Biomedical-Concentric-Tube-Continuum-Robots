function calcReachableSpace(OD, ID, k, Ls, Lc, E, modelID, nPoints, dq, simulationID)
%% This function estimates the reachable workspace for a given concentric robot
% Inputs:
%         ID: CTR Inner diameter  [m] 
%         OD: CTR Outer diameter  [m] %Replace with CTR OD, ID, k, Ls, Lc, E
%         k:  k
%         Ls: Straight section length
%         Lc: Curved section length
%         E:  Youngs Modulus
%         modelID:      identifier of the ear model 
%         nPoints:      number of points to be sampled by RRT
%         dq:           size of deltaQ, step size of rrt
%         simulationID: identifier of this simulation
fprintf('* Estimation of the reachable workspace *\n')

% add dependencies
addpath('kinematics')
addpath('utils')
addpath('utils/stlTools')
addpath('path-planning')
addpath('../anatomical-models')

%% params to keep track of
results.modelID = modelID;
results.nPoints = nPoints;
results.dq = dq;
results.simID = simulationID;

%% Part 1. Run RRT
% fprintf('Running RRT...\n')

robot = ConcentricTubeRobot(OD, ID, k, Ls, Lc, E);

% Read the configuration file to extract information about the meshes
fid = fopen(fullfile('..', 'anatomical-models', 'configurations.txt'));
text = textscan(fid, '%s %f %f %f %f %f %f %f %f %f %f %f %f');
fclose(fid);

configurations = cell2mat(text(2:end));
line_no = find(strcmp(text{1}, modelID));

path = fullfile('..', 'anatomical-models', modelID);

init_config  = configurations(line_no, 1:6);
entry_point  = configurations(line_no, 7:9);
tip_base     = configurations(line_no, 10:12);

% Calculate the transformation to the space of the CT scan
newZ = tip_base .* 1e-3 - entry_point .* 1e-3;
newZ = newZ ./ norm(newZ);
v = cross([0 0 1], newZ);
R = eye(3) + skew(v) + skew(v)^2 * (1-dot([0 0 1], newZ))/norm(v)^2;
R = R * [0 -1 0; 1 0 0; 0 0 1];
t = entry_point .* 1e-3;
T = [R t'; 0 0 0 1];

% Read the meshes from file
pathStl = fullfile('..', 'anatomical-models', modelID, 'tissue.stl');
[vertices, faces, ~, ~] = stlRead(pathStl); 
earModel.vertices = vertices*1e-3; %%Replace with our model
earModel.faces = faces; %Replace with our model
earModel.baseTransform = T; %Replace with our model

% q steps for rrt
deltaQ = ones(1,6) * dq;

% Define the robot's range of motion
%%%%%% FIGURE THIS OUT %%%%%%
minAdv = 0;
maxBend = 0.025; % [1/m] %Can be removed as kappa can be directly defined
maxKappa= 1/maxBend; %Using CalcMaxCurve.m
maxTheta= deg2rad(100); % [rad]
maxDz   = 20e-3; %[m] 
maxRot  = deg2rad(360);  % [rad]
maxAdv  = 20e-3; % [m] % dz = advancement %ARBITRARY as dependent on actuator, not in scope%
 
qBounds = [-maxKappa -maxTheta 0     0 -maxRot minAdv;
            maxKappa  maxTheta  maxDz 0  maxRot maxAdv];

robot = ConcentricTubeRobot(OD, ID, k, Ls, Lc, E)
init_config(6) = minAdv;

%% Run RRT
% To run rrt_ndim run these lines
[qListNormalized,qList,pList,aList,xList, TList, collLocs] = rrt_ndim(robot, ...
    qBounds, ...
    earModel, ... %Replace with our area model
    deltaQ, ...
    nPoints, ...
    init_config, true);

RUNTIME = toc/60;
results.runtime = RUNTIME;

fprintf('Simulation: %s completed \n', simulationID)
fprintf('RRT Runtime: %.2f minutes \t', RUNTIME)
fprintf(['Total sampled points: ' num2str(size(qList,2)) '\n']);

save([simulationID '.mat']);
end
