function makefigurevisibility(simulationID)

load([simulationID '.mat']);
fid = fopen(fullfile('..', 'anatomical-models', 'configurations.txt'));
text = textscan(fid, '%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
fclose(fid);

configurations = cell2mat(text(2:end));
line_no = find(strcmp(text{1}, modelID));

image_size   = configurations(line_no, 1:3);
voxel_size   = configurations(line_no, 4:6);


path = fullfile('..', 'anatomical-models', modelID);
load(fullfile(path, 'record.mat'));

% Read the Raw Meshes from file
pathMe = fullfile(path, 'me.mesh');
pathOs = fullfile(path, 'ossicle.mesh');
rawMeMesh = meshread(pathMe);
rawOsMesh = meshread(pathOs);

% Convert the raw meshes into objects that can be passed
% to the `patch' function
meMesh.Faces = rawMeMesh.triangles' + 1;
meMesh.Vertices = bsxfun(@times, rawMeMesh.vertices', voxel_size);
%meMesh.Vertices = meMesh.Vertices .* 1e-3;

osMesh.Faces = rawOsMesh.triangles' + 1;
osMesh.Vertices = bsxfun(@times, rawOsMesh.vertices', voxel_size);
%osMesh.Vertices = osMesh.Vertices .* 1e-3;

meMesh.FaceVertexCData = ones(size(meMesh.Vertices, 1), 1);
meMesh.LineStyle = 'none';
meMesh.FaceColor = 'flat';
meMesh.FaceAlpha = 0.4 ;
meMesh.FaceVertexCData = seenMap;
%recordnew(:);

osMesh.FaceVertexCData = ones(size(osMesh.Vertices, 1), 1);
osMesh.FaceColor = 'flat';
osMesh.FaceAlpha = 0.4 ;

figure
patch(meMesh);
xlabel('X[mm]')
ylabel('Y[mm]')
zlabel('Z[mm]')
view(-118.5, 37.74);
%legend({'Antrum', 'Epitympanum', 'Eustachian Tube', 'Hypotympanum', 'Sinus Tympanum', 'Supratubal Recess', 'Facial Recess', 'Mesotympanum'});
%title(['Simulation of endoscope with ' num2str(n) ' notches']);
set(gca,'FontSize',16);
grid on
axis equal
end