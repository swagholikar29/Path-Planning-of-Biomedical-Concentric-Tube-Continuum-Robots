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
data = readtable("Trial3_5b_3a_88deg.csv", opts);

clear opts

%% Visualize Data
animate = false;

num = length(data.Frame);
trans = zeros(4,4,num);

figure();
set(gcf, 'WindowState', 'maximized');

subplot(2,2,[1 3]);
hold on; grid on; axis equal;
title('Tip Position of Magnetic Tracker')
xlabel('X')
ylabel('Y')
zlabel('Z')
view([-90 10]);

X = data.Tx;
Y = data.Ty;
Z = data.Tz;

% h1 = scatter3(X, Y, Z, 10, data.Frame, 'Filled');

for i = 1:10:num/2
    if QuaternionsOn
        quat = table2array(data(i, 2:5));
        rot = quat2rotm(quat);
        p  = table2array(data(i, 6:8))';
    else
        rot = eul2rotm(table2array(data(i, 2:4)), "ZYX");
        p  = table2array(data(i, 5:7))';
    end
    
    T = [rot p; 0 0 0 1];
    
    triad('Matrix', T, 'scale', 1);
end

% plot of z axis over time
% subplot(2,2,2);
% hold on; grid on;
% xlabel('Time')
% ylabel('Angle (deg)')
% title('Angle from Pos')
% 
% pts = [data.Tx, data.Ty, data.Tz];

% plot of angle over time
subplot(2,2,2);
hold on; grid on;
xlabel('Time')
ylabel('Angle (deg)')
title('Angle Over Time')

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

plot(1:num, rel_angles); 

if animate
    for i = 1:5:num
        h1.XData = X(1:i);
        h1.YData = Y(1:i);
        h1.ZData = Z(1:i);
        h1.CData = data.Frame(1:i);

%         h2.XData = 1:i;
%         h2.YData = Z(1:i);
        pause(.001);
    end
end







