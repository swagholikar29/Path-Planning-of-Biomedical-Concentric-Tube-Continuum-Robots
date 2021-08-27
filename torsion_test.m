
%% CREATE PRECURVED TUBES
clear, clc, close all;
% parameters for generating tubes
ODs = [2.39e-3 1.6e-3];         % (m) outer diameters 
IDs = [2.01e-3 0];             % (m) inner diameters
E = 75e6;

precurves = [.0099e3 .0138e3];
Ls = [93.5e-3 218.5e-3];                    % (m) length of straight section
Lc = [92.3e-3 85e-3];                     % (m) length of curved section

robot = ConcentricTubeRobot(ODs, IDs, precurves, Ls, Lc, E);
robot.v = .35;

q = [10e-3   deg2rad(90);        % outermost tube
     10e-3  deg2rad(0)];     % innermost tube

robot.fwkine(q, false);

figure
degs = 90:45:270;

for i=1:length(degs)
    
    subplot(round(length(degs)/2), 2, i);
    
    robot.plotEnergyContour([deg2rad(0) deg2rad(degs(i))]);
end


% robot.plotTubes();

n = 270;
angles = zeros(n, 1);
pause(.1);

% for i = 1:n
%     
%     q = [20e-3   deg2rad(i);        % outermost tube
%          20e-3  deg2rad(0)];     % innermost tube
%     robot.fwkine(q, false);
% %     robot.animateTubes();
%     
%     arcs = robot.arcs;
%     angles(i) = rad2deg(arcs(3,2,end));
%     
% %     disp(['Input ' num2str(i) ' Output ' num2str(angles(i))])
%     pause(.01);
% end
% figure
% plot(1:n, angles)
% xlabel('Input Angle (deg)')
% ylabel('Output Angle (deg)')
% title('Angle of tubes with torision')
% 
% arcs = robot.arcs;
% disp('Curvature (1/m)');
% disp(arcs(:,1,end));
% 
% disp('Rotation (deg)');
% disp(rad2deg(arcs(:,2,end)));
