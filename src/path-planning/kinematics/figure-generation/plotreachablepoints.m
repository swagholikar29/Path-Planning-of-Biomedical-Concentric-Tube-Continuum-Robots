function plotreachablepoints(simulationID)

load([simulationID '.mat']);
figure
hold on

pathStl = fullfile('..', 'anatomical-models', modelID, 'me.stl');
[vertices, faces, ~, ~] = stlRead(pathStl);
earModel.vertices = vertices;
earModel.faces = faces;
stlPlot(earModel.vertices*1e3, earModel.faces, '');
stlPlot(osModel.vertices*1e3, osModel.faces, '', 10);

scatter3(pList(1,:)*1e3, pList(2,:)*1e3, pList(3,:)*1e3, 'filled', 'red');

axis equal
xlabel('X[mm]')
ylabel('Y[mm]')
zlabel('Z[mm]')
view(-118.5, 37.74);

%trisurf(k, pList(1,:)', pList(2,:)', pList(3,:)','FaceColor','red','FaceAlpha',0.1)

legend({'Ear Cavity', 'Ossicles', 'Reachable points'});
%title(['Reachable points with ' num2str(n) ' cutouts']);

set(gca,'FontSize',18);
end