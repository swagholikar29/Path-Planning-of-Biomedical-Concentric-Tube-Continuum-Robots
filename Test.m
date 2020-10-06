%% test

clc, clear, close all
addpath('kinematics')

%% Create Wrist

params = [];
count = 1;
x = .2:.001:.8;
for i = x
    [r, k, kj, theta, thetaj, s] = get_Curvature(i);
    params(count,:) = [k, kj, theta, thetaj, s];
    count = count + 1;
end

figure
subplot(231)
plot(x, params(:,1))
title('k')

subplot(232)
plot(x, params(:,2))
title('kj')

subplot(233)
plot(x, rad2deg(params(:,3)))
title('theta')

subplot(234)
plot(x, rad2deg(params(:,4)))
title('thetaj')

subplot(235)
plot(x, params(:,5))
title('s')

%%
[maxk, idx] = max(params(:,1));
maxk
maxg = x(idx)

[robot, newk, kj, theta, thetaj, s] = get_Curvature(maxg);
ybar = robot.ybar(1);
deltaL = (newk * robot.cutouts.h(1) * (robot.ID/2 + ybar))/(1 + newk * ybar);

q = [deltaL, 0, 0];
robot.fwkine(q, eye(4));

figure
robotModel = robot.makePhysicalModel();
X = robotModel.surface.X;
Y = robotModel.surface.Y;
Z = robotModel.surface.Z;
surf(X, Y, Z, 'FaceColor','blue');
axis equal


