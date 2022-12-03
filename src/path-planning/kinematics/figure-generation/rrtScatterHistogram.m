function rrtScatterHistogram(simulationID)
% playscatter(simulationID) creates a video animation of the endoscope and
% wrist moving to RRT points and plotting the configuration space points
% reached by the endoscope and wrist in real time (based off playresults.m)
% 
%  Authors: I. Chan <iachan@wpi.edu>
%           A. Chiluisa <ajchiluisa@wpi.edu>
%           L. Fichera  <lfichera@wpi.edu>
%  
% Last Version: 6/16/2020

load([simulationID '.mat']);

path = fullfile('..', 'anatomical-models', modelID);
pathStl = fullfile('..', 'anatomical-models', modelID, 'tissue_cropped.stl');
[vertices, faces, ~, ~] = stlRead(pathStl);
earModel.vertices = vertices*1e-3;
earModel.faces = faces;

%% Scatterplot
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,3, [1 2 4 5]);

% Visualize the robot inside the cavity
ii = 1;
numFaces = length(faces);
initialmap = ones(numFaces, 1);
h1 = stlPlot(earModel.vertices, earModel.faces, '', 0);
colormap(flipud(parula(2)));
hold on

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
scatter3(T(1,4), T(2,4), T(3,4), 100, 'k', 'filled');
hold on

robot.fwkine(qList(:,ii), T);                   % generate fwkin
robotPhysicalModel = robot.makePhysicalModel(); % generate meshes

% handler for mesh plots
h2 = surf(robotPhysicalModel.surface.Xw, ...
    robotPhysicalModel.surface.Yw, ...
    robotPhysicalModel.surface.Zw, ...
    'FaceColor','blue');
h3 = surf(robotPhysicalModel.surface.Xe, ...
    robotPhysicalModel.surface.Ye, ...
    robotPhysicalModel.surface.Ze, ...
    'FaceColor','red');
xlabel('X [m]');
ylabel('Y [m]');
zlabel('Z [m]');
title('3D Model in STL')
set(gca,'FontSize',12);

axis equal

% Endoscope Scatter Plot
subplot(2,3,3);
endo = scatter3(qList(1,1), qList(2,1), qList(3,1), 'filled');
hold on
xlabel('Curvature (Kappa) [1/m]');
ylabel('Rotation (Theta) [rad]');
zlabel('Translation (Dz) [m]');
title('Endoscope Scatterplot');
set(gca,'FontSize',12);

% Wrist Scatter Plot
subplot(2,3,6);
wrist = scatter3(qList(4,1), qList(5,1), qList(6,1), 'filled');
hold on
xlabel('Tendon Displacement (delta L) [m]');
ylabel('Rotation (Phi) [rad]');
zlabel('Advancement (Tau) [m]');
title('Wrist Scatterplot');
set(gca,'FontSize',12);

sgt = sgtitle('Sampling-based Simulation RRT Animation');
sgt.FontSize = 18;

video = VideoWriter([simulationID '.avi']); %create the video object
open(video); %open the file for writing

colorMap = zeros(numFaces, 1);
colorMap(end) = 5;
h1.FaceVertexCData = colorMap;

%for ii = 1 : size(pList, 2)
while true
    robot.fwkine(qList(:,ii), T);
    robotPhysicalModel = robot.makePhysicalModel();
    
    % Visibility
    v = visibleMap(:,ii);
    colorMap = ones(numFaces, 1);
    colorMap(logical(v)) = 0;
    h1.FaceVertexCData = colorMap;
    
    % Physical Model
    h2.XData = robotPhysicalModel.surface.Xw;
    h2.YData = robotPhysicalModel.surface.Yw;
    h2.ZData = robotPhysicalModel.surface.Zw;
    
    h3.XData = robotPhysicalModel.surface.Xe;
    h3.YData = robotPhysicalModel.surface.Ye;
    h3.ZData = robotPhysicalModel.surface.Ze;
    title(['Pose ' num2str(ii) ' of ' num2str(size(pList, 2))]);
    
    % Endoscope (Kappa, Theta, L)
    endo.XData = qList(1,1:ii);
    endo.YData = qList(2,1:ii);
    endo.ZData = qList(3,1:ii);
    endo.CData = [repmat([0 0.4470 0.7410], ii-1, 1); 1 0 0];
    endo.SizeData = [36 * ones(1,ii-1), 100];
    
    % Wrist (Tau, Phi, Delta L)
    wrist.XData = qList(4,1:ii);
    wrist.YData = qList(5,1:ii);
    wrist.ZData = qList(6,1:ii);
    wrist.CData = [repmat([0 0.4470 0.7410], ii-1, 1); 1 0 0];
    wrist.SizeData = [36 * ones(1,ii-1), 100];
    
    
    ii = ii + 1;
    if ii > size(pList, 2), break; end
    
    currFrame = getframe(gcf);
    writeVideo(video,currFrame);
  
end
close(video);
%% Histograms (comment and uncomment as needed)
% figure
% subplot(2, 3, 1)
% histogram(qListNormalized(1,:)) 
% title('Histogram of Endoscope Curvature')
% xlabel('Kappa [1/m]')
% 
% subplot(2, 3, 2)
% histogram(qListNormalized(3,:))
% title('Histogram of Endoscope Translation')
% xlabel('s [m]')
% 
% subplot(2, 3, 3)
% histogram(qListNormalized(2,:))
% title('Histogram of Endoscope Base Rotation')
% xlabel('Theta [rad]')
% 
% subplot(2, 3, 4)
% histogram(qListNormalized(4,:))
% title('Histogram of Wrist Tendon Displacement')
% xlabel('Delta l [m]')
% 
% subplot(2, 3, 5)
% histogram(qListNormalized(5,:))
% title('Histogram of Wrist Base Rotation')
% xlabel('Phi [rad]')
% 
% subplot(2, 3, 6)
% histogram(qListNormalized(6,:))
% title('Histogram of Wrist Advancement')
% xlabel('Tau [m]')

%% Scatterplots
% dim = 7;
% figure
% scatter(qListNormalized(1,:), qListNormalized(2,:), 'b.')
% title('Scatterplot of RRT Points (1 vs. 2)')
% xlabel('Dimension 1')
% ylabel('Dimension 2')
% axis equal
% 
% figure
% scatter(qListNormalized(1,:), qListNormalized(dim,:), 'r.')
% title(['Scatterplot of RRT Points (1 vs. ', num2str(dim), ')'])
% xlabel('Dimension 1')
% ylabel(['Dimension ', num2str(dim)])
% axis equal
end