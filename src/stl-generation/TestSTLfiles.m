% Plot the STL files to manually determine the entry point (X, Y, Z) in
% the Larynx.

% Select the larynx model you want to visualize acording with the folder
% name.

clear all; clc;

modelID = 'larynx8b';

%
p = [ 0.8       0.8       1.0;    % green
      0.2422    0.1504    0.6603; % purple
      0.9769    0.9839    0.0805]; % yellow
% Files that we are going to read
file1 = 'tissue.stl';
file2 = 'tissue_cropped.stl';

% Obtain and read the files from the directory
Stl1 = fullfile('..', '..', 'anatomical-models', modelID, file1);
Stl2 = fullfile('..', '..', 'anatomical-models', modelID, file2);
[V1,F1,Name1] = stlRead(Stl1);
[V2,F2,Name2] = stlRead(Stl2);


%% Plot the STL files 

% % Plot the larynx for calcReachableSpace function 
% figure
% stlPlot(V1,F1,Name1,p(2,:));
% title(modelID)
% %title(upper(file1(1:6)));

% Plot the larynx model for calcVisibleArea function
figure
stlPlot(V2,F2,Name2,p(2,:));
title(modelID)
%title('TISSUE CROPPED');