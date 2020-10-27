%% Image Analyzer Script
% goes through series of images and applies the filter to get the tube
% then uses circle fit on the tube points

close all; clear all; clc; 
addpath('images');

%% Tube Parameters
OD = 3;     % (mm) outer diameter
Lc = 50;    % (mm) curved section length

%% Fix Image
img_path = 'Nylon/Nylon_Tube_3mm/T5_3mm.jpg';
% img_path = 'test.jpg';
img = imread(img_path);        % load image
img = imrotate(img, 90);
[bw, rgb] = maskBlackTubes(img);    % mask out not black
bw2 = bwareaopen(bw, 200000);
bw3 = imfill(bw2, 'hole');
bw3 = imfill(bw3, 'hole');

% Plot masked images
if false
    figure
    subplot(131)
    imshow(bw)
    subplot(132)
    imshow(bw2)
    subplot(133)
    imshow(bw3)
end

%% Get bottom edge of tube

[yi, xi] = find(bw3);   % indicies of white points 
yi = -yi;               % flip y 

ye = [];
xe = [];
delay = 1;
count = 1;

% iterate through all possible indices of X
for i = min(xi) + delay:max(xi)
    % skip if i is not in x
    if ~ismember(xi, i)
        continue;
    end
    
    % check if arc starts to curve back up
    newy = min(yi(xi == i));
    if i > 1000 && length(ye) > delay+1
        % check curr y against previous ones with delay
        if abs(newy) < abs(ye(count - delay))
%             break;
        end
    end
    
    % store min y for inner side of tube
    xe(count) = i;
    ye(count) = min(yi(xi == i));
    count = count + 1;
end




%% ------ Calc R by first calculating bending angle
% ** Assumes known length of bending section
[xc,yc,Re] = circle_fit(xe, ye);
center = [xc, yc];      % coords for center of circle
[inner] = calcBendingAngle([xe(1) ye(1)], [xe(end) ye(end)], center);
unscaled_R = (inner/Lc) * 1e3;

%% Calc R 
[xc,yc,Re_px] = circle_fit(xe,ye);
center = [xc, yc];      % coords for center of circle

% Calc bending angle of curve
[inner_ang start_ang ustart utip] = calcBendingAngle([xe(1) ye(1)], [xe(end) ye(end)], center);

% calculate scale in mm/px
scale = (Lc/inner_ang - OD/2)/Re_px;

% scaled parameters
xe = xe*scale;
ye = ye*scale;
Re = Re_px * scale;
R = Re + OD/2;
center = center * scale;

%% Plot results
% create circle 
% th = linspace(pi/2 - start_ang,pi/2 - inner_ang - start_ang,100)';
th = linspace(0, 2*pi, 100);
xfit = Re*cos(th)+center(1);
yfit = Re*sin(th)+center(2);

x_real = R*cos(th)+center(1); 
y_real = R*sin(th)+center(2);

figure
hold on;
scatter(xi*scale,yi*scale, 'b', 'DisplayName', 'Mask of Tube')
scatter(xe,ye,'r', 'DisplayName', 'Fit Edge of Tube')
plot(x_real,y_real,'-.', 'LineWidth', 1.5, 'DisplayName', 'Curve Fit of Tube');
plot(xfit,yfit,'-.', 'LineWidth', 2, 'DisplayName', 'Curve Fit of Edge');
quiver([center(1);center(1)],[center(2); center(2)], [ustart(1); utip(1)] * 100, [ustart(2); utip(2)] * 100,...
    'DisplayName', 'Bending Angle');

title([img_path 'Measured fitted and true circles'])
% text('Location', 'southwest',yc,sprintf('R=%g mm \nCurvatue=%g 1/m \nAngle=%g rad',R,1000/R, inner_ang))
txt = sprintf('R= %.2fmm \nCurvatue= %.2f 1/m \nAngle= %.2fdeg',R,1000/R, rad2deg(inner_ang));
annotation('textbox', [.6 0 .5 .3], 'String', txt, 'FitBoxToText', 'on')
legend('Location', 'southwest')
xlabel('X (mm)')
ylabel('Y (mm)')
axis equal; grid on;

disp(txt);
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



