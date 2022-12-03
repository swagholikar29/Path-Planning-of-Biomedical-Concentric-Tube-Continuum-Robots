function makeVisibilityFig(simulationID)
% Produces figures were the visible area of a model is highlighted
% 
% simulationID = identifier of this simulation
%
% Author: Jesse F. d'Almeida  <jfdalmeida@wpi.edu>

load([simulationID '.mat']);
fid = fopen(fullfile('..', 'anatomical-models', 'configurations.txt'));
text = textscan(fid, '%s %f %f %f %f %f %f %f %f %f %f %f %f');
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
meMesh.faces = faces;
meMesh.vertices = vertices .* 1e-3;

plotReachable = true;

if exist('visibleMapTotal', 'var')
    v = logical(visibleMapTotal);
else
    v = zeros(length(faces),1);
end

if plotReachable
    figure('Name', simulationID);
    subplot(121);
    stlPlot(meMesh.vertices * 1e3, meMesh.faces, 'Visibility', v);
%     hold on, axis equal
%     scatter3(collLocs(:,1)*1e3, collLocs(:,2)*1e3, collLocs(:,3)*1e3, 'filled', 'green');
    
    subplot(122);
    stlPlot(meMesh.vertices * 1e3, meMesh.faces, 'Visibility', v);
    hold on, axis equal
    scatter3(pList(1,:)*1e3, pList(2,:)*1e3, pList(3,:)*1e3, 'filled', 'red');
else
    figure
    v = zeros(length(faces),1);
    stlPlot(meMesh.vertices * 1e3, meMesh.faces, 'Visibility', v);
    hold on, axis equal
    scatter3(pList(1,:)*1e3, pList(2,:)*1e3, pList(3,:)*1e3, 'filled', 'red');
end
end