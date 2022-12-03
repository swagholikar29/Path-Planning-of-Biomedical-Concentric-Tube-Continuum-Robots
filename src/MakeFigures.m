clear; clc; close all;

simulationID = 'larynx1-dq-0.06-100pts';
% 
load([simulationID '.mat']);
% fid = fopen(fullfile('..', 'anatomical-models', 'configurations.txt'));
% text = textscan(fid, '%s %f %f %f %f %f %f %f %f %f %f %f %f');
% fclose(fid);
% 
% configurations = cell2mat(text(2:end));
% line_no = find(strcmp(text{1}, modelID));
% 
% path = fullfile('..', 'anatomical-models', modelID);

% % Read the Raw Meshes from file
% pathMe = fullfile(path, 'tissue_cropped.stl');
% [vertices, faces, ~, ~] = stlRead(pathMe);
% numverts = size(vertices, 2);
% 
% % Convert the raw meshes into objects that can be passed
% % to the `patch' function
% meMesh.faces = faces;
% meMesh.vertices = vertices .* 1e-3;
% meMesh.vertices1 = vertices .*1e-3;
% 
% 
if exist('visibleMapTotal', 'var')
    v = logical(visibleMapTotal);
else
    v = zeros(length(faces),1);
end

figure('Name', simulationID);
stlPlot(meMesh.vertices * 1e3, meMesh.faces, 'Visibility', v);
axis equal 

if LarynxModel == "L1"
    view(-45, 35) 
    xlim([17 40])
    zlim([-58 -20])
elseif LarynxModel == "L2"
    view(-135,35)
    xlim([-57 -20])
    ylim([-42 -17])
    zlim([-60 -20])    
end

set(gca, 'FontName', 'CMU Serif', 'fontsize', 20);
title('Steerable Fiber');

