function h = plotAllTubes(tubes, q)
%PLOTALLTUBES plots meshes of all tubes without deformation
%   Must run fwkine on all models before running this plotter
%   INPUT
%       tubes: [N] array of Precurved objects
%   OUTPUT
%       h: handles for plots
numTubes = length(tubes);
h = zeros(1,numTubes);
colors = distinguishable_colors(numTubes);

%% Plot Entire Model
% figure('Name', 'Individual Tubes');
% hold on

for i = 1:numTubes
    figure
    tube = tubes(i);
    p = q(i, 1);
    rot = q(i,2);
    
    arc = [0 rot p;
           tube.precurve 0 tube.Lc];
    tube.fwkine(arc);
    
    model = tube.makePhysicalModel();
    
    % mesh model of the tube
    h(i) = surf(model.surface.X, model.surface.Y, model.surface.Z,...
        'FaceColor', colors(i,:));
    
    % create a triad (coord frame) for the transformations
    trans = tubes(i).transformations;
    for jj = 1:size(trans,3)
        triad('Matrix', trans(:,:,jj), 'scale', 5e-3);
    end

    axis('image');
    view([135 30]);
    grid on;

    camlight('headlight');
    material('dull');
    axis equal
%     xlabel('X (m)');
%     ylabel('Y (m)');
%     zlabel('Z (m)');
    title('Individual Tubes No Deformation');
end
end

