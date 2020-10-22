%% Image Analyzer Script
% goes through series of images and applies the filter to get the tube
% then uses circle fit on the tube points

close all; clear all; clc; 
addpath('images');

%% Tube Parameters
OD = 4;     % (mm) outer diameter
Lc = 50;    % (mm) curved section length

%% Fit circle
img = imread('cropped.jpg');        % load image
[bw, rgb] = maskBlackTubes(img);    % mask out not black
bw2 = imfill(bw, 'hole');
bw3 = bwareaopen(bw2, 10);

% indicies of white points 
[yi, xi] = find(bw3);
yi = -yi;

minYs = [];
minXs = [];
lasty = 0;

% iterate through all possible indices of X
for i = min(xi):2500
    % skip if i is not in x
    if ~ismember(xi, i)
        continue;
    end
    
    % check if arc starts to curve back up
    newy = min(yi(xi == i));
%     if newy > lasty
%         continue;
%     end
    
    lasty = newy;
    
    % store min y for inner side of tube
    minXs(i) = i;
    minYs(i) = min(yi(xi == i));
end
% set the scale based on tube size
%   set by comparing width of img at base to known tube diameter
x1 = min(xi);           % initial x's
y1 = yi(xi == x1 );     % all initial y's
ydiff = abs(max(y1) - min(y1));
scale = OD/ydiff;       % scale in mm/px

x = minXs;
y = minYs;

[xc,yc,Re] = circle_fit(x,y);
center = [xc, yc];      % coords for center of circle
[inner] = calcBendingAngle([minXs(1) minYs(1)], [minXs(end) minYs(end)], center);
unscaled_R = (inner/Lc) * 1e3

x = x*scale;
y = y*scale;

[xc,yc,Re] = circle_fit(x,y);
center = [xc, yc];      % coords for center of circle

R = Re + OD/2;          % Radius of curvature at center is R at edge offset by tube radius

% Calc bending angle of curve
[inner_ang start_ang ustart utip] = calcBendingAngle([x(1) y(1)], [x(end) y(end)], center);

%% Plot results
% create circle 
th = linspace(pi/2 - start_ang,pi/2 - inner_ang - start_ang,100)';
xfit = Re*cos(th)+xc;
yfit = Re*sin(th)+yc;

x_real = R*cos(th)+xc; 
y_real = R*sin(th)+yc;

figure
hold on;
scatter(xi*scale,yi*scale, 'DisplayName', 'Mask of Tube')
scatter(x,y,'r', 'DisplayName', 'Fit Edge of Tube')
plot(x_real,y_real,'-.', 'LineWidth', 1.5, 'DisplayName', 'Curve Fit of Tube');
plot(xfit,yfit,'-.', 'LineWidth', 2, 'DisplayName', 'Curve Fit of Edge');
quiver([center(1);center(1)],[center(2); center(2)], [ustart(1); utip(1)] * 10, [ustart(2); utip(2)] * 10,...
    'DisplayName', 'Bending Angle');

title('Measured fitted and true circles')
% text('Location', 'southwest',yc,sprintf('R=%g mm \nCurvatue=%g 1/m \nAngle=%g rad',R,1000/R, inner_ang))
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



