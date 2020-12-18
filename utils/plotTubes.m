function h = plotTubes(tubes)
%PLOTTUBES plots meshes of all tubes in the array
%   Must run fwkine on all models before running this plotter
%   INPUT
%       tubes: array of Precurved objects
%   OUTPUT
%       h: handles for plots

plotBackbone = false;
plotIndividualDeformations = false;


numTubes = length(tubes);
h = zeros(1,numTubes);
colors = distinguishable_colors(numTubes);

%% Plot Entire Model
figure('Name', 'Precurved Tubes');
hold on

for i = 1:numTubes
    model = tubes(i).makePhysicalModel();
    
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
%     zlim([-.01 .08]);
%     ylim([-.05 .05]);
%     xlim([-.05 .05]);
    xlabel('X (m)');
    ylabel('Y (m)');
    zlabel('Z (m)');
    title('Concentric Precurved Tubes with Deformation');
end

%% on separate figures plot the backbone points of each tube 
% (for debugging)
if plotBackbone
    for i = 1:numTubes
        figure
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

%% Plot individually deformed tubes
if plotIndividualDeformations
figure('Name', 'Individual Deformed Tubes');
    for i = 1:numTubes
        subplot(1,numTubes,i);
        
        model = tubes(i).makePhysicalModel();
    
        % mesh model of the tube
        h(i) = surf(model.surface.X, model.surface.Y, model.surface.Z,...
            'FaceColor', colors(i,:));

        axis('image');
        view([-135 35]);
        grid on;

        camlight('headlight');
        material('dull');
        axis equal
        xlabel('X (m)');
        ylabel('Y (m)');
        zlabel('Z (m)');
        title(['Deformed Tube ' num2str(i)]);
    end
end

