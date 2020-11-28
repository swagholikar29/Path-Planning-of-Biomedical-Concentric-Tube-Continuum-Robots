%% Image Analyzer Script
% goes through series of images and applies the filter to get the tube
% then uses circle fit on the tube points

close all; clear all; clc; 
addpath('images');


%% Make image name from params

% Tube Parameters
Lc = 50;    % (mm) curved section length
OD = [6 4];     % (mm) outer diameter
material = 'PA12';
modifier = 'c';
trial = 0;

switch material
    case 'PA12'
        short = 'PA';
    case 'Nylon'
        short = 'NY';
    case 'Durable'
        short = 'DU';
end
f = fullfile('C:/', 'Users', 'jfd42', 'OneDrive', 'Documents', 'WPI-DESKTOP-AGKBP1C', 'Year 5- Masters', 'Thesis', 'trials\');

if length(OD) == 1
    img_path = sprintf('%s\\Tube_%dmm\\%s%dmm_T%d%s.JPG', material, OD, short, OD, trial, modifier)
else
    img_path = sprintf('%s\\Tube_%dmm_%dmm\\%s%dmm_%dmm_T%d%s.JPG', material, OD(1), OD(2), short, OD(1), OD(2), trial, modifier)
    OD = OD(1);
end

path = [f img_path];

% list of ODs that do not work well with full circ fit and instead need to fit edge
ODforEdgeFit = [2 3 4 6]; 
fitEdge = ismember(OD, ODforEdgeFit);   % boolean for switch algs

%% Fix Image
img = imread(path);        % load image
img = imrotate(img, 0);
[bw, rgb] = maskBlackTubes(img);    % mask out not black
bw2 = bwareaopen(bw, 200000);
bw3 = imfill(bw2, 'hole');
bw3 = imfill(bw3, 'hole');

plotmask = false;
plotcircs = true;

% Plot masked images
if plotmask
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

if fitEdge
    disp("Finding Edge and fitting circle");
    [bot_edge, top_edge] = getEdge(xi, yi);
    [txc,tyc,tRe_px] = circle_fit(top_edge(1,:), top_edge(2,:));
    [xc,yc,Re_px] = circle_fit(bot_edge(1,:), bot_edge(2,:));
    xe = bot_edge(1,:);
    ye = bot_edge(2,:);
else
    disp("Fitting Circle");
    [bot_edge, top_edge] = getEdge(xi, yi);
    [xc,yc,Re_px] = circle_fit(xi, yi);
    xe = xi;
    ye = yi;
end

%% Calc R 
center = [xc, yc];      % coords for center of circle

% Calc bending angle of curve
[inner_ang start_ang ustart utip] = calcBendingAngle([xe(1) ye(1)], [xe(end) ye(end)], center);

scale = calcScaleByOD(OD, bot_edge, top_edge)

% calculate scale in mm/px
% scale = (Lc/inner_ang - OD/2)/Re_px;

% scaled parameters
xe = xe*scale;
ye = ye*scale;
Re = Re_px * scale;

center = center * scale;

% if the fit is on the edge then need to get actual radius by adding radius
% of tube
if fitEdge
    tRe = tRe_px * scale;
    R = Re + OD/2;
else
    R = Re;
end

%% Plot results
% create circle 
th = linspace(0, 2*pi, 100);
xfit = Re*cos(th)+center(1);
yfit = Re*sin(th)+center(2);

x_real = R*cos(th)+center(1); 
y_real = R*sin(th)+center(2);



txt = sprintf('R= %.2fmm \nCurvatue= %.2f 1/m \nAngle= %.2fdeg',R,1000/R, rad2deg(inner_ang));
disp(txt);

if plotcircs
    figure
    hold on;
    scatter(xi*scale,yi*scale, 'b', 'DisplayName', 'Mask of Tube')
    scatter(xe,ye,'r', 'DisplayName', 'Fit Edge of Tube')
    plot(x_real,y_real,'-.', 'LineWidth', 2, 'DisplayName', 'Curve Fit of Tube');
    
    if fitEdge
        x_top = tRe*cos(th)+txc*scale; 
        y_top = tRe*sin(th)+tyc*scale;
        scatter(top_edge(1,:)*scale, top_edge(2,:)*scale, 'DisplayName', 'Top Edge of Tube')
        plot(xfit,yfit,'-.', 'LineWidth', 1.5, 'DisplayName', 'Curve Fit of Edge');
        plot(x_top,y_top,'-.', 'LineWidth', 1.5, 'DisplayName', 'Top Curve');
    end

%     quiver([center(1);center(1)],[center(2); center(2)], [ustart(1); utip(1)] * 100, [ustart(2); utip(2)] * 100,...
%         'DisplayName', 'Bending Angle');
    
    title([img_path])
    % text('Location', 'southwest',yc,sprintf('R=%g mm \nCurvatue=%g 1/m \nAngle=%g rad',R,1000/R, inner_ang))
    
    annotation('textbox', [.6 0 .5 .3], 'String', txt, 'FitBoxToText', 'on')
    legend('Location', 'southwest')
    xlabel('X (mm)')
    ylabel('Y (mm)')
    axis equal; grid on;

end

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

function [bot_edge, top_edge] = getEdge(xi, yi)
    % CIRCFITEDGE finds the points along the bottom edge of the tube
    %INPUTS:
    %   xi = x indices of the masked tube
    %   yi = y indices of the masked tube
    %OUTPUTS:
    %   bot_edge = [2xN] coords of the bottom tube edge
    %   top_edge = [2xN] coords of the top tube edge
    
    bot_edge = [];
    top_edge = [];
    delay = 10;
    startdelay = 100;
    count = 1;

    % iterate through all possible indices of X
    for i = min(xi) + startdelay:max(xi)
        % skip if i is not in x
        if ~ismember(xi, i)
            continue;
        end

        % check if arc starts to curve back up
        newy = min(yi(xi == i));
        if i > 2500 && size(bot_edge, 2) > delay+1
            % check curr y against previous ones with delay
            if abs(newy) < abs(bot_edge(2, count - delay))
                break;
            end
        end

        % store min y for inner side of tube
        bot_edge(1, count) = i;
        bot_edge(2, count) = min(yi(xi == i));
        
        top_edge(1, count) = i;
        top_edge(2, count) = max(yi(xi == i));
        
        count = count + 1;
    end
end


function scale = calcScaleByOD(OD, bot_edge, top_edge)
    % CALCSCALEBYOD finds the mm/px scale of the image by finding the
    % diameter, in pixels, of the tube and using the known OD
    %INPUTS:
    %   OD = (mm) outer diameter of tube
    %   bot_edge = [2xN] coords of the bottom tube edge
    %   top_edge = [2xM] coords of the top tube edge
    %OUTPUTS:
    %   scale = (mm/px) scale of tube
    
    n = 20;  % number of points to sample
    min_dists = [];
%     topRange = 200;
    pts = floor(linspace(1, size(bot_edge, 2), n + 1));  % find points along tube to sample
    pts = pts(2:end-1);
    
    % at each point on bottom edge, calc min dist to points along top edge
    for i = 1:length(pts)
        bot = bot_edge(:,pts(i));
        top = top_edge(:,pts(i));
        curr_min_dist = norm(bot - top);     % init with point in same column
        
        % go through points in top edge to find min dist
        for j = 1:size(top_edge,2)
            d = norm(bot - top_edge(:,j));
            if d < curr_min_dist
                curr_min_dist = d;
            end
            
        end
        min_dists(i) = curr_min_dist;
    end
    
    avg_dist = mean(min_dists);
    scale = OD/avg_dist;
    
end

