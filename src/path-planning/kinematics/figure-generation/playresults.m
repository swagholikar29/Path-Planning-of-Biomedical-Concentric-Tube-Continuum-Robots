function playresults(simulationID)

load([simulationID '.mat']);

path = fullfile('..', 'anatomical-models', modelID);
pathStl = fullfile(path, 'me.stl');
[vertices, faces, ~, ~] = stlRead(pathStl);
earModel.vertices = vertices;
earModel.faces = faces;


figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,2,1);

% Visualize the robot inside the cavity
ii = 1;
h1 = stlPlot(earModel.vertices, earModel.faces, '');
stlPlot(osModel.vertices, osModel.faces, '');
hold on

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

robot.fwkine(qList(:,ii), T);
robotPhysicalModel = robot.makePhysicalModel();
h2 = surf(robotPhysicalModel.surface.X, ...
    robotPhysicalModel.surface.Y, ...
    robotPhysicalModel.surface.Z, ...
    'FaceColor','blue');
xlabel('X [m]');
ylabel('Y [m]');
zlabel('Z [m]');
set(gca,'FontSize',14);

axis equal

subplot(1,2,2);
h3 = scatter3(qList(1,1), qList(2,1), qList(3,1), 'filled');
hold on
xlabel('Tendon Displacement [m]');
ylabel('Rotation [rad]');
zlabel('Advancement [m]');
set(gca,'FontSize',14);


sgt = sgtitle('Sampling-based Simulation of Middle Ear Endoscopy');
sgt.FontSize = 18;

video = VideoWriter([simulationID '.avi']); %create the video object
open(video); %open the file for writing

%for ii = 1 : size(pList, 2)
while true
    robot.fwkine(qList(:,ii), T);
    robotPhysicalModel = robot.makePhysicalModel();
    
    h2.XData = robotPhysicalModel.surface.X;
    h2.YData = robotPhysicalModel.surface.Y;
    h2.ZData = robotPhysicalModel.surface.Z;
    title(['Pose ' num2str(ii) ' of ' num2str(size(pList, 2))]);
    
    h3.XData = qList(1,1:ii);
    h3.YData = qList(2,1:ii);
    h3.ZData = qList(3,1:ii);
    h3.CData = [repmat([0 0.4470 0.7410], ii-1, 1); 1 0 0];
    h3.SizeData = [36 * ones(1,ii-1), 100];
    
    %drawnow
    
    %     fprintf('Press "n" to move forward or "p" to move back.\n')
    %     fprintf('Press any other key to stop testing and generate the reachable workspace.\n\n')
    %
    %     while ~waitforbuttonpress, end
    %     k = get(gcf, 'CurrentCharacter');
    
    %     switch k
    %         case 'p'
    %             ii = ii - 1;
    %             if ii < 1, ii = 1; end
    %         case 'n'
    ii = ii + 1;
    if ii > size(pList, 2), break; end
    %if ii > 3000, break; end
    %         otherwise
    %             break
    %     end
    currFrame = getframe(gcf);
    writeVideo(video,currFrame);
    %ii
end

close(video);
end