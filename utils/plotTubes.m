function h = plotTubes(tubes)
%PLOTTUBES plots meshes of all tubes in the array
%   Must run fwkine on all models before running this plotter
%   INPUT
%       tubes: array of Precurved objects
%   OUTPUT
%       h: handles for plots

numTubes = length(tubes);
h = zeros(1,numTubes);
colors = distinguishable_colors(numTubes);

figure('Name', 'Precurved Tubes');
hold on

for i = 1:numTubes
    model = tubes(i).makePhysicalModel();
    
    h(i) = surf(model.surface.X, model.surface.Y, model.surface.Z,...
        'FaceColor', colors(i,:));

    axis('image');
    view([-135 35]);
    grid on;

    camlight('headlight');
    material('dull');
    axis equal
%     zlim([-.01 .08]);
%     ylim([-.05 .05]);
%     xlim([-.05 .05]);
    xlabel('X (m)');
    ylabel('Y (m)');
    zlabel('Z (m)');
end
end

