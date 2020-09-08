%% Script to create series of precurved concentric tubes and evaluate their kinematics and mechanics

close all, clear, clc
addpath('kinematics');

tube = Precurved(5e-3, 90, 30e-3, 30e-3);

q = [0 30e-3];
tube.fwkine(q);
model = tube.makePhysicalModel();

figure
surf(model.surface.X, model.surface.Y, model.surface.Z)

axis('image');
view([-135 35]);
grid on;

camlight('headlight');
material('dull');
axis equal
zlim([-.01 .08]);
ylim([-.05 .05]);
xlim([-.05 .05]);
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
