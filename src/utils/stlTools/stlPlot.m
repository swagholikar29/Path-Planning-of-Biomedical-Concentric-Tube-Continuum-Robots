function h = stlPlot(v, f, name, seenMap, trueColors)
%STLPLOT is an easy way to plot an STL object
% Outputs:
%       h: patch object handle
% Inputs: 
%       V: (Nx3) array of vertices
%       F: (Mx3) array of faces
%       NAME: (string) name of the object, that will be displayed as a title
%       seenMap: (Mx1) array of the visibility of faces (accepts either
%           logical or double values)
%       trueColors: (Px3) array mapping of value to a rgb color based on
%           index

%figure, hold on
object.vertices = v;
object.faces = f;

if nargin < 4
    h = patch(object,'FaceColor', [0.8 0.8 1.0], ...
        'EdgeColor',       'none', 'EdgeALpha', .3,       ...
        'FaceLighting',    'gouraud',     ...
        'FaceAlpha', 0.4, ...
        'AmbientStrength', 0.15);
else
    paint = true;
    if ~exist('trueColors', 'var')
        paint = false;
    end
    
    colors = zeros(length(f), 1);
    seenMap = double(seenMap);  % incase seenmap is just logical, convert 
    
    if ~paint
        colors = seenMap;
    else
        % use input colormap for face colors
        colors = trueColors(seenMap+1, :);
    end
    h = patch(object,'FaceVertexCData', colors, ...
        'FaceColor', 'flat', ...
        'EdgeColor',       'none',        ...
        'FaceLighting',    'gouraud',     ...
        'FaceAlpha', 0.6, ...
        'AmbientStrength', 0.15);    
end

xlabel('X [mm]'), ylabel('Y [mm]'), zlabel('Z [mm]');

% Add a camera light, and tone down the specular highlighting
camlight('headlight');
material('dull');
zlim([-.05 .1]);
ylim([-.1 .1]);
xlim([-.1 .1]);

% Fix the axes scaling, and set a nice view angle
axis('image');
view([-135 35]);
grid on;
title(name);
