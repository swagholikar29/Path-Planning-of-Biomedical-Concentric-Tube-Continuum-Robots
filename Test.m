%% test

clc, clear, close all
addpath('kinematics')

%% Find best depth of cut

% creates wrist at each depth of cut to find max curve
params = [];

counti = 1;
countj = 1;
x = .001: .0001: .03;
lenX = length(x);

rs = cell(lenX);
for i = x
    [r, k, kj, theta, thetaj, s] = get_Curvature(i);
    params(counti, :) = [k, kj, theta, thetaj, s];
    rs{counti} = r;

    counti = counti + 1;
end

figure
subplot(121)
plot(x, params(:,1))
grid on;
title('Curvature of the Wrist')
xlabel('max strain')
ylabel('curvature (1/m)');

%% Plot wrist at max
[maxk, idx] = max(params(:,1));
maxBendingRadius = 1/maxk * 1e3;
sAtMax = params(idx,5);  % max s
robot = rs{idx};        % best robot

% [robot, newk, kj, theta, thetaj, s] = get_Curvature(maxParam);
ybar = robot.ybar(1);
h = robot.cutouts.h(1); 
deltaL = h - 2 * (1 / maxk - robot.ID/2) * sin((maxk * h) / (2 * (1 + ybar * maxk)));

disp(['Min bending radius: ' num2str(maxBendingRadius)  ' mm ']);
disp(['Max curvature ' num2str(maxk)  ' 1/m']);
disp(['Length of bending section: ' num2str(sAtMax*1e3)  ' mm']);
disp(['Best param: ' num2str(x(idx))  ' mm']);
q = [deltaL, 0, 10e-3];
robot.fwkine(q, eye(4));

% Precurved(robot.OD, robot.ID, maxk, 10e-3, s)

subplot(122)
robotModel = robot.makePhysicalModel();
X = robotModel.surface.X;
Y = robotModel.surface.Y;
Z = robotModel.surface.Z;
surf(X, Y, Z, 'FaceColor','blue');
title('Model of Wrist at Max Curvature')
axis equal


