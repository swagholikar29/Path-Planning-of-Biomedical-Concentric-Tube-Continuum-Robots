%% Figure Maker
clear all; close all; clc;

%% Tubes Plot

% parameters for generating tubes
% ODs = [5.6 3.9].*1e-3;         % (m) outer diameters 
ODs = [4 3 2] .* 1e-3;
IDs = ODs - 1.2e-3;             % (m) inner diameters
E = 1.515e9;

precurves = [15 30 45];
Ls = [75e-3 140e-3 140e-3];                    % (m) length of straight section
Lc = [50e-3 50e-3 50e-3];                     % (m) length of curved section
robot = ConcentricTubeRobot(ODs, IDs, precurves, Ls, Lc, E);

q = [20e-3   deg2rad(0);        % outermost tube
     50e-3  deg2rad(180);
     80e-3      deg2rad(0)];     % innermost tube
% q = [20e-3 deg2rad(-90)];

robot.fwkine(q, false);
figure
hold on
robot.plotTubes();
view([180 0]);

set(gca,'FontName','Times New Roman','FontSize',14);
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
set(gca,'ztick',[])
set(gca,'zticklabel',[])

figure
precurves = [0 0 0];
Ls = [75e-3 140e-3 140e-3];                    % (m) length of straight section
Lc = [50e-3 50e-3 50e-3];                     % (m) length of curved section
robot = ConcentricTubeRobot(ODs, IDs, precurves, Ls, Lc, E);
q = [5e-3   deg2rad(0);        % outermost tube
     30e-3  deg2rad(0);
     60e-3      deg2rad(0)];     % innermost tube

robot.fwkine(q, false);
robot.plotTubes();

view([180 0]);
set(gca,'FontName','Times New Roman','FontSize',14);
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
set(gca,'ztick',[])
set(gca,'zticklabel',[])

%% Box and Whisker Plots
inplane_diffs = [1.4	0.03	-1.01	-0.02;
1.35	0.56	1.65	-0.82;
0.26	-0.09	1.59	-0.25;
3.14	0.7	1.67	0.22;
1.05	0.33	2.25	0.19;
2.04	-1.08	2.51	0.19;
1.51	-0.15	2.26	-0.13;
1.73	-1.54	2.68	-0.14;
0.26	0.55	2.12	-0.26;
-1.73	0.89	1.96	0.17];
ip_labels = {'Trial A', 'Trial B', 'Trial C', 'Trial D'};
% boxplot(inplane_diffs,  'Labels', ip_labels);

trans_diff = [-9.564273243	3.470903631	0.593855906;
nan	-1.167228556	1.364334649;
-8.999975766	-0.73930991	0.622033008;
-3.725831482	-6.511664682	-1.260069097;
6.913280217	-1.811130588	2.340032896;
nan	-5.076911099	0.86992403];
trans_labels = {'Link 1', 'Link 2', 'Link 3'};
boxplot(trans_diff, 'Labels', trans_labels);
			
% Create ylabel
ylabel('Difference in Bending Radius (mm)','FontSize',16,'FontName','Times New Roman');

box('on');
grid('on');
% Set the remaining axes properties
set(gca,'FontName','Times New Roman','FontSize',14);

%% Error Bar Plots in-plane bending
figure
avgs = [1.499090909	0.720909091	1.858181818	0.249090909];
std_dev = [0.777413597	0.595428868	0.575433721	0.200655949];
stderr = std_dev ./ sqrt(11);
x = 1:length(avgs);

bar(x,avgs);
hold on; grid on; 
box('on');
er = errorbar(x, avgs, stderr, stderr);
title("Error in Radius for In-Plane Bending");
er.Color = [.5 0 0];       
er.LineStyle = 'none';
er.LineWidth = 2;

ylabel('Difference in Bending Radius (mm)','FontSize',16,'FontName','Times New Roman');
xlabel('Trial','FontSize',16,'FontName','Times New Roman');
set(gca,'FontName','Times New Roman','FontSize',14);

%% Error Bar Plot for translation
figure
avgs = [9.412798474	3.129524744	1.175041598];
std_dev = [4.18341725	2.107511176	0.596791183];
stderr = std_dev ./ sqrt(6);
x = 1:length(avgs);

bar(x,avgs);
hold on; grid on; 
box('on');
er = errorbar(x, avgs, stderr, stderr);
title("Error in Radius for Translation");
er.Color = [.5 0 0];       
er.LineStyle = 'none';
er.LineWidth = 2;

ylabel('Difference in Bending Radius (mm)','FontSize',16,'FontName','Times New Roman');
xlabel('Link','FontSize',16,'FontName','Times New Roman');
set(gca,'FontName','Times New Roman','FontSize',14);

%% Error bar Rotation

figure
avgs = [11.3982 2.319 2.2936 0.5184];
std_dev = [0.351142421	0.203 0.351142421	0.203];
stderr = std_dev ./ sqrt(11);
x = 1:length(avgs);

bar([1 2], avgs(1:2), 'FaceColor', 'red');
hold on
bar([3 4], avgs(3:4), 'FaceColor', 'blue');
hold on; grid on; 
box('on');
er = errorbar(x, avgs, stderr, stderr);
title("Error in Rotation Angle");
er.Color = [.5 0 0];       
er.LineStyle = 'none';
er.LineWidth = 2;


ylabel('Difference in Angle (deg)','FontSize',16,'FontName','Times New Roman');
xlabel('Combination','FontSize',16,'FontName','Times New Roman');
set(gca,'FontName','Times New Roman','FontSize',14);

%% Error bar full kin

figure
avgs = [0.432234552	7.745616792	3.204756934	3.612248036];
std_dev = [0.015612724	0.642538199	1.572836814	0.277872568];
stderr = std_dev ./ sqrt(length(avgs));
x = 1:length(avgs);

bar(x, avgs);
hold on
% bar([3 4], avgs(3:4), 'FaceColor', 'blue');
hold on; grid on; 
box('on');
er = errorbar(x, avgs, stderr, stderr);
title("Error in Tip Position");
er.Color = [.5 0 0];       
er.LineStyle = 'none';
er.LineWidth = 2;


ylabel('Difference in Tip Position (mm)','FontSize',16,'FontName','Times New Roman');
xlabel('Configuration','FontSize',16,'FontName','Times New Roman');
set(gca,'FontName','Times New Roman','FontSize',14);

%% Nitinol Stress Strain Curve

epsl = 0.018;
epsp = 0.085;

Elin = 40e9;
Ese = 0.08*Elin;
Ep = 1e9;
beta = 0.01;

sigl = Elin * epsl;
sigp = sigl + (epsp - epsl)*Ese;

n = 100;
e1 = linspace(0, epsl, n);
e2 = linspace(epsl, epsp, n);
e3 = linspace(epsp, 0.14, n);

y1 = e1 * Elin;
y2 = sigl + (e2 - epsl)*Ese;
y3 = sigp + Ep*exp(-beta./(e3-epsp));

figure
plot([e1 e2 e3], [y1 y2 y3]./1e6, 'LineWidth',2,'Color',[1 0 0])
grid on;

% Create ylabel
ylabel('Stress (MPa)','FontSize',16,'FontName','Times New Roman');

% Create xlabel
xlabel('Strain','FontSize',16,'FontName','Times New Roman');

box('on');
grid('on');
% Set the remaining axes properties
set(gca,'FontName','Times New Roman','FontSize',14);

%% Stress Strain Curve
m = readmatrix('CSVs/dogbone.csv');

n = 5;
stress = m(1:end-n, 4);
strain = m(1:end-n, 5);


plot(strain, stress, 'LineWidth',2,'Color',[1 0 0])
grid on;

% Create ylabel
ylabel('Stress (MPa)','FontSize',16,'FontName','Times New Roman');

% Create xlabel
xlabel('Strain (%)','FontSize',16,'FontName','Times New Roman');
title('Experimental Stress Stain Curve of PA12 Black');
box('on');
grid('on');
% Set the remaining axes properties
set(gca,'FontName','Times New Roman','FontSize',14);

%% Deformation 

x1 = [0 10 20 50 75 100];

pa12_4 = [0
3.965788135
7.670867062
9.015484246
11.98888813
10.08386813];


nylon_4 = [0
4.895222483
-1.180217006
2.568204633
8.129553437
8.083302728];

durable_4 = [0
8.837491993
12.2163155
15.60825289
21.87101775
24.58005096];

x2 = [0 10 20 30 40 50 60 70 80 90 100];

pa12_3 = [0.00
5.44
5.02
6.51
5.21
6.10
6.42
5.18
6.42
6.43
6.07];

nylon_3 = [0.00
1.51
1.21
2.55
1.67
3.18
1.92
2.34
2.89
2.32
1.22];

durable_3 = [0.00
5.56
6.99
10.55
11.00
9.47
14.95
16.77
16.42
20.57
21.21];


labels4 = {'PA12', 'Nylon', 'Durable'};
figure
subplot(2, 1, 1)
plot(x1, pa12_4, x1, nylon_4, x1, durable_4, 'LineWidth',2)
legend(labels4);
grid on;
axis equal;
ylim([-5, 35]);
xlim([0 100]);
% Create ylabel
ylabel('Deformed Radius (mm)','FontSize',16,'FontName','Times New Roman');

% Create xlabel
xlabel('Number of Cycles','FontSize',16,'FontName','Times New Roman');
title('Plastic Deformation of 4 mm Diameter Tubes');

box('on');
grid('on');
% Set the remaining axes properties
set(gca,'FontName','Times New Roman','FontSize',14);

subplot(2,1,2);
plot(x2, pa12_3, x2, nylon_3, x2, durable_3, 'LineWidth',2)
legend(labels4);
grid on;
axis equal;
ylim([-5, 35]);
xlim([0 100]);
% Create ylabel
ylabel('Deformed Radius (mm)','FontSize',16,'FontName','Times New Roman');

% Create xlabel
xlabel('Number of Cycles','FontSize',16,'FontName','Times New Roman');
title('Plastic Deformation of 3 mm Diameter Tubes');

box('on');
grid('on');
% Set the remaining axes properties
set(gca,'FontName','Times New Roman','FontSize',14);







