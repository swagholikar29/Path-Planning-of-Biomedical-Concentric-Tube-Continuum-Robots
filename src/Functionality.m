clear; close; clc
load('larynx1-dq-0.06-10pts.mat');

path = fullfile('..', 'anatomical-models', modelID);
pathStl = fullfile('..', 'anatomical-models', modelID, 'tissue_cropped.stl');
[vertices, faces, ~, ~] = stlRead(pathStl);
earModel.vertices = vertices*1e-3;
earModel.faces = faces;

%% Scatterplot
figure('units','normalized','outerposition',[0 0 1 1])
% subplot(2,3, [1 2 4 5]); -- ADD LATER WITH THE REST OF THE REMOVED LINES

% Visualize the robot inside the cavity
ii = 1;
numFaces = length(faces);
colorMap = zeros(numFaces, 1);
% 
% % color indices for painting the larynx
% colors = [0.2422    0.1504    0.6603;   % purple
%           0.9769    0.9839    0.0805];   % yellow


%color indices for painting the larynx
colors = [1 0.5 0.45;   % skin pink
          0.9769    0.9839    0.0805];   % yellow

h1 = stlPlot(earModel.vertices, earModel.faces, '', colorMap, colors);
hold on

% plot the visibility rays coming out of the tip
plotQuiv = false;

if plotQuiv
    axis manual
    vp = quiv(ii).vp;
    rays = quiv(ii).dispR;
    q1 = quiver3(vp(1,:), vp(2,:), vp(3,:), rays(1,:), rays(2,:), rays(3,:));
    yl = ylim;
    xl = xlim; 
    zl = zlim;
end

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
    'FaceColor','#808080');
h3 = surf(robotPhysicalModel.surface.Xe, ...
    robotPhysicalModel.surface.Ye, ...
    robotPhysicalModel.surface.Ze, ...
    'FaceColor','black');
xlabel('X [m]');
ylabel('Y [m]');
zlabel('Z [m]');
title('3D Model in STL')
set(gca,'FontSize',12);

% this sets the limits of the axis. When plotting quiver, if not hard
% defined, it likes to change the limits based on the plots so the figure
% moves around a bit
if plotQuiv
    ylim(yl);
    xlim(xl);
    zlim(zl);
    axis manual
end

v = logical(visibleMap(:,1));

axis off
zlim([-38e-3 -20e-3])
view(0,0)

video = VideoWriter([simulationID '.avi']); %create the video object
video.FrameRate = 10;
open(video); %open the file for writing

while true
    
    % Quiver of rays
    if plotQuiv
        % scale down the rays to be in units of m
        vp = quiv(ii).vp * 1e-3;
        rays = quiv(ii).dispR * 1e-3;
        
        %plot rays 
        q1.XData = vp(1,:);
        q1.YData = vp(2,:);
        q1.ZData = vp(3,:);
        q1.UData = rays(1,:);
        q1.VData = rays(2,:);
        q1.WData = rays(3,:);
    end
    
    robot.fwkine(qList(:,ii), T);
    robotPhysicalModel = robot.makePhysicalModel();
    
%     % Visibility
%     if exist('visibleMap', 'var')
%         v = any([v visibleMap(:,ii)], 2); % calculate the row-wise logical OR
%         colorMap = zeros(numFaces, 3);
%         colorMap = colors(logical(v)+1, :);
%         h1.FaceVertexCData = colorMap;
%     end    
    
    % Physical Model
    h2.XData = robotPhysicalModel.surface.Xw;
    h2.YData = robotPhysicalModel.surface.Yw;
    h2.ZData = robotPhysicalModel.surface.Zw;
    
    h3.XData = robotPhysicalModel.surface.Xe;
    h3.YData = robotPhysicalModel.surface.Ye;
    h3.ZData = robotPhysicalModel.surface.Ze;
%    title('Laser Fiber Configuration Space');
%    sgt.String = ['Sampling-based Exploration of the Larynx - Tissue Coverage Estimation - Pose ' num2str(ii) ' of ' num2str(nPoints)];
    
    % Endoscope (Kappa, Theta, L)
    endo.XData = qList(1,1:ii);
    endo.YData = qList(2,1:ii);
    endo.ZData = qList(3,1:ii);
%     endo.CData = [repmat([0 0.4470 0.7410], ii-1, 1); 1 0 0];
    endo.SizeData = [36 * ones(1,ii-1), 100];
    
    % Wrist (Tau, Phi, Delta L)
    wrist.XData = qList(4,1:ii);
    wrist.YData = qList(5,1:ii);
    wrist.ZData = qList(6,1:ii);
%     wrist.CData = [repmat([0 0.4470 0.7410], ii-1, 1); 1 0 0];
    wrist.SizeData = [36 * ones(1,ii-1), 100];
    
    savefig(['figures/IROS/' 'sim-' num2str(ii) '.fig']);
    saveas(gcf,['figures/IROS/' 'sim-' num2str(ii), '.svg']);    
    saveas(gcf,['figures/IROS/' 'sim-' num2str(ii), '.jpg']); 
    
    ii = ii + 1;
    if ii > size(pList, 2), break; end
    
    currFrame = getframe(gcf);
    writeVideo(video,currFrame);
  
end
close(video);