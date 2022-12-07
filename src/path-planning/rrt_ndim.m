 function [qListNormalized,qList,pList,aList,xList, TList, collLocs] = rrt_ndim(robot, qbounds, model, deltaQ, nPoints, init_config)
% RRT implements the basic Rapidly-Exploring Random Trees algorithm for a
% generic continuum robot
%
%   robot: (robot obj) instance of Robot object = ConcentricTubeRobot
%   qbounds: (2x6) bounds for each parameter
%       row 1: min bounds
%       row 2: max bounds
%   model: (model obj) STL model of region
%   deltaQ: (1x6) step size of each parameter
%   nPoints: (int) amount of points to explore
%   init_config: (1x6) initial configuration of robot

if nargin < 3
    collisionDetection = false;
    nPoints = 1000;
else
    collisionDetection = true;
end

% collisionDetection = false;

n = length(deltaQ);     % number of parameters

% function to scale the normalized configurations
minBounds = qbounds(1,:)';
maxBounds = qbounds(2,:)';
scaleByBounds = @(q) q.*(maxBounds - minBounds) + minBounds;

% normalize init_config
norm_init_config = (init_config.' - minBounds) ./ (maxBounds - minBounds);

% take care of any anomalies
check_nan = isnan(norm_init_config);
for i = 1:size(norm_init_config)
    if check_nan(i) == 1 || abs(norm_init_config(i)) > 1
        norm_init_config(i) = 0;
    end
end

% initialize the tree and the starting point
qListNormalized = zeros(n, nPoints);
qListNormalized(:,1) = norm_init_config;
qList = zeros(n, nPoints);
qList(:,1) = init_config;
pList = zeros(3, nPoints);
aList = zeros(3, nPoints);
xList = zeros(3, nPoints);   % x-direction for each points
TList = zeros(4,4, nPoints); % tip transformations for each point

T_robot_in_env = model.baseTransform;

robot.fwkine(qList(:,1), T_robot_in_env);   % configure fwkin
T = robot.transformations(:,:,end);         % transformation to tip
pList(:,1) = T(1:3,4,end);                  % initial pose
aList(:,1) = T(1:3,3,end);                  % initial approach
xList(:,1) = T(1:3,1,end);                  % initial x direction

% iteratively build the tree
hw = waitbar(0, 'Sampling the configuration space. Please wait...');

% determine the bounds for the larynx model
minX = min(model.vertices(:,1));
maxX = max(model.vertices(:,1));
minY = min(model.vertices(:,2));
maxY = max(model.vertices(:,2));
minZ = min(model.vertices(:,3));
maxZ = max(model.vertices(:,3));
boundsCount = 0;

% vertices of faces in the model
v1 = model.vertices(model.faces(:,1),:);
v2 = model.vertices(model.faces(:,2),:);
v3 = model.vertices(model.faces(:,3),:);

jj = 1;
lastTime = 0;
dispInterval = 60*5; % time in seconds between displaying loop info
collLocs = [];
nColl = 0;
tic
while true
    % Generate a random point, identify the closest point in the tree
    % and move towards the new point
    now = toc;
    if (now - lastTime) >= dispInterval
        lastTime = now;
        estTime =((nPoints/jj) * now - now)/60;
        fprintf('On loop %d of %d \t %.2f minutes into simulation with about %.2f minutes left\n',...
            jj, nPoints, now/60, estTime)
    end
    
    qRand = rand(n,1);      % random point of n-dim
    qNearest = nearestVertex(qRand, qListNormalized, jj);
    qNewNorm = move(qNearest, qRand, deltaQ);
    
    % Scale up the point
    qNew = scaleByBounds(qNewNorm);
    
    robot.fwkine(qNew, T_robot_in_env);     % fwkin of new configuration
    T = robot.transformations(:,:,end);
    
    % Check for potential collisions
    if collisionDetection
        robotPM = robot.makePhysicalModel();        % generate meshes
        
        % list of points in mesh of wrist and endoscope
        testpts = [robotPM.surface.X robotPM.surface.Y robotPM.surface.Z;
                           robot.pose(:,end)'];
        
        collision = intriangulation(model.vertices, ...
            model.faces, testpts);
        
        collision = sum(collision);
        
        if collision > 0
            
            % cannot operate if initial configuration collides
            if jj == 1
                warning('Initial configuration results in collision')
                return
            end
            
            % display the collision number only if colliding a lot
            nColl = nColl + 1;
            if mod(nColl, 200) == 0 && nColl ~=0
%                 fprintf('Collisions detected at this configuration: %d \n', nColl);
            end          
            continue;
        end
        
       % check if within the model space
       if Ptip(1) < minX || Ptip(1) > maxX || ...
               Ptip(2) < minY || Ptip(2) > maxY || ...
               Ptip(3) < minZ || Ptip(3) > maxZ
          boundsCount = boundsCount + 1;
          continue;
       end
    end
    
    % If no collision, add this point to the tree
    qListNormalized(:,jj) = qNewNorm;
    
    qList(:,jj) = qNew;
    
    pList(:,jj) = T(1:3,4,end);
    aList(:,jj) = T(1:3,3,end);
    xList(:,jj) = T(1:3,1,end);
    TList(:,:,jj) = T(:,:,end);
    
    jj = jj + 1;
    
    nColl = 0;
    
    if jj > nPoints
        break;
    end
    
    waitbar(jj/nPoints, hw, 'Sampling the configuration space. Please wait...');
end

close(hw);

qListNormalized = qListNormalized(:,1:jj-1);
qList = qList(:,1:jj-1);
pList = pList(:,1:jj-1);
aList = aList(:,1:jj-1);
xList = xList(:,1:jj-1);
TList = TList(:,:,1:jj-1);
end