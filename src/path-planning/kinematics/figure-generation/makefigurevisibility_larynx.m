function makefigurevisibility(simulationID)

load([simulationID '.mat']);
fid = fopen(fullfile('..', 'anatomical-models', 'configurations.txt'));
text = textscan(fid, '%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
fclose(fid);

configurations = cell2mat(text(2:end));
line_no = find(strcmp(text{1}, modelID));

path = fullfile('..', 'anatomical-models', modelID);

% Read the Raw Meshes from file
pathMe = fullfile(path, 'tissue_cropped.stl');
[vertices, faces, ~, ~] = stlRead(pathMe);
numverts = size(vertices, 2);

% Convert the raw meshes into objects that can be passed
% to the `patch' function
meMesh.Faces = faces;
meMesh.Vertices = vertices .* 1e-3;

plotReachable = true;

if plotReachable
    figure('Name', simulationID)
    subplot(121)
    stlPlot(meMesh.Vertices * 1e3, meMesh.Faces, 'Visibility', visibleMapTotal);
    
    subplot(122)
    stlPlot(meMesh.Vertices * 1e3, meMesh.Faces, 'Visibility', visibleMapTotal);
    hold on, axis equal
    scatter3(pList(1,:)*1e3, pList(2,:)*1e3, pList(3,:)*1e3, 'filled', 'red');
    set(gca,'FontSize',16);
else
    figure
    stlPlot(meMesh.Vertices * 1e3, meMesh.Faces, 'Visibility', visibleMapTotal);
end
end