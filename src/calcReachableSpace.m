function calcReachableSpace(u,h,w,ID,OD,modelID, nPoints, dq, useWrist, simulationID) %Replace with(OD, ID, k, Ls, Lc, E, modelID, nPoints, dq, simulationID)
%% This function estimates the reachable workspace for a given endoscope configuration
% Inputs:
%         u: length of uncut sections   [m] %Remove from arguments
%         h: length of cut sections     [m] %Remove from arguments
%         w: width  of the cut sections [m] %Remove from arguments
%         ID: endoscope inner diameter  [m] %Replace with CTR OD, ID, k, Ls, Lc, E
%         OD: endoscope outer diameter  [m] %Replace with CTR OD, ID, k, Ls, Lc, E
%
%         modelID:      identifier of the ear model %%%% replace with our model %%%%
%         nPoints:      number of points to be sampled by RRT
%         dq:           size of deltaQ, step size of rrt
%         useWrist:     (boolean) false- bounds changed to no wrist
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

% Define the endoscope model
n = length(u); % number of notches
cutouts.w = w;
cutouts.u = u; % converting to meters
cutouts.h = h; % converting to meters
cutouts.alpha = zeros(1,n);  

robot = EndoWrist(ID, OD, n, cutouts); %Change robot || robot = ConcentricTubeRobot(OD, ID, k, Ls, Lc, E)

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

%%%Remove%%%
% now slide the endoscope back by its length, so that all the different
% designs start exploring from the same point
wrist_len_offset = sum(cutouts.u) + sum(cutouts.h); % len of wrist bend section
endo_len_offset = 13.4e-3 + 28.2e-3;                % tip len + bend section len


Tz = eye(4);
Tz(3,4) = -wrist_len_offset - endo_len_offset;
Tz(3,4) = - endo_len_offset;
T = T * Tz;
%%%Remove%%%

% Read the meshes from file
pathStl = fullfile('..', 'anatomical-models', modelID, 'tissue.stl');
[vertices, faces, ~, ~] = stlRead(pathStl); 
earModel.vertices = vertices*1e-3; %%Replace with our model
earModel.faces = faces; %Replace with our model
earModel.baseTransform = T; %Replace with our model

% q steps for rrt
%        [kappa   endoBaseRot dz tendonDisp wristBaseRot wristAdvancement]
deltaQ = ones(1,6) * dq;

% Define the robot's range of motion
%%%%%% FIGURE THIS OUT %%%%%%
minAdv = -sum(cutouts.h) - sum(cutouts.u) + 3e-3; %Re-define to zero
maxBend = 0.025; % [1/m] %Can be removed as kappa can be directly defined
maxKappa= 1/maxBend; %Replace with calcMaxCurve.m
maxTheta= deg2rad(100); % [rad]
maxDz   = 20e-3; %[m]   %ARBITRARY as dependent on actuator, not in scope%
maxDisp = sum(cutouts.h); % [m]  %ARBITRARY as dependent on actuator, not in scope%
maxRot  = deg2rad(360);  % [rad]
maxAdv  = 20e-3; % [m] % dz = advancement %ARBITRARY as dependent on actuator, not in scope%

% Normal Bounds TO BE IGNORED
if useWrist
    qBounds = [-maxKappa, -maxTheta, 0,     0,      -maxRot, minAdv;
                maxKappa,  maxTheta, maxDz, maxDisp, maxRot, maxAdv];
else 
% No wrist bending Bounds %TO BE USED% %as we have no wrist
    qBounds = [-maxKappa -maxTheta 0     0 -maxRot minAdv;
               maxKappa  maxTheta  maxDz 0  maxRot maxAdv];
end

robot = EndoWrist(ID, OD, n, cutouts); %Replace with ctr || robot = ConcentricTubeRobot(OD, ID, k, Ls, Lc, E)

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
