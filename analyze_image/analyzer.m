%% Image Analyzer Script
% goes through series of images and applies the filter to get the tube
% then uses circle fit on the tube points

close all; clear all; clc; 
addpath('images');

%% Tube Parameters
OD = 4;     % (mm) outer diameter
Lc = 50;    % (mm) curved section length

%% Read Image and Adjust
img = imread('cropped.jpg');        % load image
[bw, rgb] = maskBlackTubes(img);    % mask out not black
bw2 = bwareaopen(bw, 1000);
bw3 = imfill(bw2, 'hole');           % 
% figure
% subplot(131)
% imshow(bw)
% subplot(132)
% imshow(bw2)
% subplot(133)
% imshow(bw3)

%% Get edge of the circle
[y_bw, x__bw] = find(bw3);  % indicies of white points 
y_bw = -y_bw;               % flip image for some reason

Y = [];
X = [];
delay = 0;     % to catch when bottom edge starts to curve up around the tip

% iterate through all possible indices of X
% for i = min(x__bw)+delay:max(x__bw)
for i = 10:2500
    % skip if i is not in x
    if ~ismember(x__bw, i)
        continue;
    end
    
    % check if arc starts to curve back up and stop there
    newy = min(y_bw(x__bw == i));
%     if newy == 0
%         continue;
%     end
%     if length(Y) > delay+1 && i > 1000
%         % check curr y against previous ones with delay
%         if abs(newy) < abs(Y(i-delay))
%             break;
%         end
%     end
    
    % store min y for inner side of tube
    X(i) = i;
    Y(i) = min(y_bw(x__bw == i));
end

%% ------ Calc R by first calculating bending angle
% ** Assumes known length of bending section
% [xc,yc,Re] = circle_fit(X,Y);
% center = [xc, yc];      % coords for center of circle
% [inner] = calcBendingAngle([X(1) Y(1)], [X(end) Y(end)], center);
% unscaled_R = (inner/Lc) * 1e3


%% ----- Calc R by scaling down
% Assumes known OD and scaling factor is correct
%   set the scale based on tube size by
%   comparing width of img at base to known tube diameter

x1 = min(x__bw);           % initial x's
y1 = y_bw(x__bw == x1 );     % all initial y's
ydiff = abs(max(y1) - min(y1));
scale = OD/ydiff;       % scale in mm/px

% scale x and y
X = X*scale;
Y = Y*scale;

[xc,yc,Re] = circle_fit(X,Y);
center = [xc, yc];      % coords for center of circle
R = Re + OD/2;          % Radius of curvature at center is R at edge offset by tube radius

% Calc bending angle of curve
[inner_ang start_ang ustart utip] = calcBendingAngle([X(1) Y(1)], [X(end) Y(end)], center);

%% Plot results
% create circle 
% th = linspace(pi/2 - start_ang,pi/2 - inner_ang - start_ang,100)';
th = linspace(0, 2*pi, 100);
xfit = Re*cos(th)+xc;
yfit = Re*sin(th)+yc;

x_real = R*cos(th)+xc; 
y_real = R*sin(th)+yc;

% plots
figure
hold on;
% scatter(x__bw*scale,y_bw*scale, 'DisplayName', 'Mask of Tube')
scatter(X,Y,'r', 'DisplayName', 'Fit Edge of Tube')
plot(x_real,y_real,'-.', 'LineWidth', 1.5, 'DisplayName', 'Curve Fit of Tube');
plot(xfit,yfit,'-.', 'LineWidth', 2, 'DisplayName', 'Curve Fit of Edge');
quiver([center(1);center(1)],[center(2); center(2)], [ustart(1); utip(1)] * 10, [ustart(2); utip(2)] * 10,...
    'DisplayName', 'Bending Angle');

title('Measured fitted and true circles')
txt = sprintf('R= %.2fmm \nCurvatue= %.2f 1/m \nAngle= %.2fdeg',R,1000/R, rad2deg(inner_ang));
annotation('textbox', [.6 0 .5 .3], 'String', txt, 'FitBoxToText', 'on')
legend('Location', 'southwest')
xlabel('X (mm)')
ylabel('Y (mm)')
axis equal; grid on;


%% Helper Functions
function [inner_ang start_ang ustart utip] = calcBendingAngle(firstpt, lastpt, center)
    % get vectors
    start = [firstpt(1), firstpt(2)] - center;   % vec of first point on edge
    tip = [lastpt(1), lastpt(2)] - center; % vec of last point on edge

    % calc angle between
    ustart = start/norm(start);
    utip = tip/norm(tip);
    inner_ang = acos(dot(ustart, utip));

    % ang from yaxis to start
    start_ang = acos(dot(ustart, [0 1]));
end



