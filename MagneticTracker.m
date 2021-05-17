%% Convert CVS from Aurora Magnetic Tracker into usable rotations and visualize

close all, clear, clc;

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 7);

QuaternionsOn = false;

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
if QuaternionsOn
    opts.VariableNames = ["Frame", "Q0", "Qx", "Qy", "Qz", "Tx", "Ty", "Tz"];
else
    opts.VariableNames = ["Frame", "Rz", "Ry", "Rx", "Tx", "Ty", "Tz"];
end

opts.VariableTypes = repmat("double", 1, length(opts.VariableNames));

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
% data = readtable("CSVs/fullkin/fullkinK.csv", opts);
data = readtable("CSVs/Trial2_5a_3b_90deg.csv", opts);

clear opts

%% Plot Transforms
animate = false;

num = length(data.Frame);
trans = zeros(4,4,num);

figure();
set(gcf, 'WindowState', 'maximized');

% subplot(2,2,[1 3]);
% hold on; grid on; axis equal;
% title('Tip Position of Magnetic Tracker')
% xlabel('X (mm)')
% ylabel('Y (mm)')
% zlabel('Z (mm)')
% view([-90 10]);
% set(gca,'fontsize', 18);
% 
% X = data.Tx;
% Y = data.Ty;
% Z = data.Tz;

% h1 = scatter3(X, Y, Z, 10, data.Frame, 'Filled');

% for i = 1:10:num/2
%     if QuaternionsOn
%         quat = table2array(data(i, 2:5));
%         rot = quat2rotm(quat);
%         p  = table2array(data(i, 6:8))';
%     else
%         rot = eul2rotm(table2array(data(i, 2:4)), "ZYX");
%         p  = table2array(data(i, 5:7))';
%     end
%     
%     T = [rot p; 0 0 0 1];
%     
%     triad('Matrix', T, 'scale', 1);
% end

% plot of z axis over time
% subplot(2,2,2);
% hold on; grid on;
% xlabel('Time')
% ylabel('Angle (deg)')
% title('Angle from Pos')
% 
% 

%% plot x y z positions
t = (1:num)/100;
dt = (1:num-1)/100;

% get relative position data
pts = [data.Tx, data.Ty, data.Tz, data.Rx];
initial_Ts = mean(pts(1:30, :));
pts = pts - initial_Ts;

% find flat regions of x (stopped parts)
dx = diff(pts(:,4));
dx0_idx = find(abs(dx) <= .0015);

% subplot(1,2,1);
hold on;
grid on;

nCols = size(pts,2);
for i=1:nCols
    plot(t, pts(:,i));
end

for i=1:nCols
    scatter(dx0_idx/100, pts(dx0_idx, i), 10, 'filled');
end
labels = {'X data', 'Y data', 'Z data', 'Angle', 'X flat', 'Y flat', 'Z flat', 'Angle flat'};
% labels = {'X data', 'Y data', 'Z data', 'X flat', 'Y flat', 'Z flat'};
legend(labels);

% find mean of each region
% go through idx
last = dx0_idx(1);
flat = [last];
avgs = zeros(0,nCols);
for i=2:length(dx0_idx)
    
    curr = dx0_idx(i);
    if curr - last <= 20            % if idx is next to last add to flat region
        flat(end+1) = curr;
    elseif length(flat) < 3       % check if theres a large flat region
        flat = [];
    else                           % find mean if at end of good flat region
        avgs(end+1, :) = mean(pts(flat, :));
        flat = [];
    end
    last = curr;
end
disp(labels);
disp(avgs);

%% plot of angle over time
% subplot(1,2,2);

t = (1:num)/100;
figure
hold on; grid on;
xlabel('Time (s)')
ylabel('Theta (deg)')
title('Angle Over Time')
set(gca,'fontsize', 18);

if QuaternionsOn
    quat0 = table2array(data(10, 2:5));
    rot0 = quat2rotm(quat0);
    z0 = rot0(1:3, 3);

    start = 200;
    fin = num-200;
    quat = table2array(data(start:fin, 2:5));
    rot = quat2rotm(quat);
    z = rot(1:3, 3, :);

    z0s = repmat(z0, 1, 1, length(z));
    
    angles = rad2deg(acos(dot(z, z0s)));
    angles = reshape(angles, length(angles), 1);
    
else
    angles = data.Rx;
    
end

angles(angles < 0) = 360 + angles(angles < 0);
initial_angle = mean(angles(1:30));

rel_angles = angles - initial_angle;

plot(t, rel_angles); 
% 
% if animate
%     for i = 1:5:num
%         h1.XData = X(1:i);
%         h1.YData = Y(1:i);
%         h1.ZData = Z(1:i);
%         h1.CData = data.Frame(1:i);
% 
% %         h2.XData = 1:i;
% %         h2.YData = Z(1:i);
%         pause(.001);
%     end
% end







