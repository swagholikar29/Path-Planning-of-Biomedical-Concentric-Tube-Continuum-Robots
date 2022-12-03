function calcVisibleArea(simulationID, alg, laserOffsetAngle)
%% This function runs a ray casting procedure to estimate the visible surface attained by a 
%  translaryngeal endoscope inside the larynx
%  
%  IMPORTANT: Before running this function, you need to run
%  `calcReachableSpace' and then pass this function the same
%  simulationID
%
%  Inputs:
%       simulationID: identifier of this simulation
%       alg - ('hpr' or 'mcrc') to select between which algorithm is used (Hidden
%       Point Removal or Monte Carlo Ray Casting)
%       laserOffsetAngle - [deg] angle from z that laser projects at
%
% Authors: Loris Fichera 
%          Jesse F. d'Almeida  <jfdalmeida@wpi.edu>
% 
% Last Revision: 7/1/2020

%% Load in Simulation
fprintf('* Estimation of the visible surface. Using %s algorithm*\n', alg)

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

% Convert the raw meshes into objects that can be passed
% to the `patch' function
meMesh.faces = faces;
meMesh.vertices = vertices .* 1e-3;
% meMesh.Vertices = meMesh.Vertices ;

meMesh.FaceVertexCData = ones(size(meMesh.vertices, 1), 1);
meMesh.LineStyle = 'none';
meMesh.FaceColor = 'flat';
meMesh.FaceAlpha = 0.4 ;

% Calculate the visual range- list of faces [3xM] that are visible
visibleFaces = [];

% empty quiver for when no faces are visible
q1.vp = [0;0;0];
q1.dispR = [0;0;0];

% empty outputs to fill
quiv = repmat(q1, 1, nPoints);          
visibleMap = zeros(length(faces),1);

hw = waitbar(0, 'Calculating the visibility map. Please wait...');

%% RUN VISIBILITY: for each point, get list of visibile faces

% use angle offset from z-axis to calculate the approach vector by
% transforming to the tip
laserAng = deg2rad(laserOffsetAngle);           % convert laser into degrees
laserVec = [sin(laserAng) 0 cos(laserAng) 0]';  % unit vector of laser direction (0 at end so that it can be transformed)

% for each point, calculate visibility 
for jj = 1 : size(pList,2)
    
    % get current transformation to tip
    T_tip = TList(:,:,jj);
    approach4 = T_tip * laserVec;                   % approach vector with extra element from transform
    approach = approach4(1:3);                      % split 3D vector

    [m, q] = visibilitymap(pList(:,jj), approach, meMesh, alg);
    visibleMap(:,jj) = m;   % visibility of faces for this point
    
    % add if there are quivers, if not add empty so no rays will be
    % displayed
    if q.vp
        quiv(jj) = q;
    else
        quiv(jj) = q1;
    end
    
    waitbar(jj/nPoints, hw, 'Calculating the visibility map. Please wait...');
end
close(hw);

% Sum rows of map to get amount for each face
numFaces = length(faces);
visibleMapTotal = sum(visibleMap, 2);            % sum rows of visible map
visFACES = sum(logical(visibleMapTotal));        % total visible face
percVISFACES = (visFACES/numFaces) * 1e2;

%% GATHER RETULS: calculate area visible
visible_area = seenArea(meMesh, visibleMapTotal);
allFaces = ones(numFaces,1);
total_area = seenArea(meMesh, allFaces);
percVISAREA = visible_area/total_area * 100;
visAREA = visible_area * 1e6;

fprintf('Faces Visible: %d \tPercent of Total Faces: %.2f%% \n',...
    visFACES, percVISFACES)
fprintf('Visible Surface Area: %.2f mm^2 \tPercent of Total Surface Area: %.2f%%\n\n',...
    visAREA, percVISAREA);

%% OUTPUTS
results.visFaces = visFACES;
results.percFaces = percVISFACES;
results.visArea = visAREA;
results.percArea = percVISAREA;

save([simulationID '.mat'], '-v7.3');


end