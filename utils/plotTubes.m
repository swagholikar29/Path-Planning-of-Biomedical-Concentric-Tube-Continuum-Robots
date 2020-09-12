function h = plotTubes(tubes)
%PLOTTUBES plots meshes of all tubes in the array
%   Must run fwkine on all models before running this plotter
%   INPUT
%       tubes: array of Precurved objects
%   OUTPUT
%       h: handles for plots

plotBackbone = 0;

numTubes = length(tubes);
h = zeros(1,numTubes);
colors = distinguishable_colors(numTubes);

figure('Name', 'Precurved Tubes');
% subplot(2,3, [1 2 4 5]);
hold on

for i = 1:numTubes
    model = tubes(i).makePhysicalModel();
    
    h(i) = surf(model.surface.X, model.surface.Y, model.surface.Z,...
        'FaceColor', colors(i,:));
    
    trans = tubes(i).transformations;
    for jj = 1:size(trans,3)
        triad('Matrix', trans(:,:,jj), 'scale', 5e-3);
    end

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

if plotBackbone
    for i = 1:numTubes
        figure
        % subplot(2,3,3);
        hold on;
        b = tubes(i).robotModel.backbone;
        p = tubes(i).pose;
        hs(i) = scatter3(b(1,:), b(2,:), b(3,:), 'filled', 'MarkerFaceColor', colors(i,:));
        scatter3(p(1,:), p(2,:), p(3,:), 'filled');
        trans = tubes(i).transformations;
        for jj = 1:size(trans,3)
            triad('Matrix', trans(:,:,jj), 'scale', 5e-3);
        end


        grid on;
        axis equal;
    end
end
end

