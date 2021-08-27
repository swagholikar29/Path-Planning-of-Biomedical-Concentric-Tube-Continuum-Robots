%% Script to create series of precurved concentric tubes and evaluate their kinematics and mechanics
close all, clear, clc
addpath('kinematics');

%% CREATE PRECURVED TUBES
% get resultant curvatures from list of tubes
Rs = [67.36796562	40.61657128];


Ks = 1 ./ (Rs * 1e-3);         % (1/m) curvatures
straight = zeros(length(Ks),1);

% Ks = [straight Ks];

ODs = [5.6e-3 3.9e-3];         % (m) outer diameters 

numRs = size(Ks,1);
nTubes = size(Ks,2);
% resultant_k = zeros(size(Ks,1),1);
resultant_k = [];
resultant_ang = [];

for i = 1:numRs
    % parameters for generating tubes
    IDs = ODs - 1.2e-3;

    Ls = 100-3;                    % (m) length of straight section
    Lc = 50e-3;                     % (m) length of curved section

    % make the precurved tubes
    for j = 1:nTubes
        tubes(j) = Precurved(ODs(j), IDs(j), Ks(i,j), Ls, Lc, 1.9e9);
    end

    %% Calculate tube deformations and forward kinematics

    % -------------change these parameters-------------
    % Joint Parameters of the tube
    % translation  rotation
    q = [0e-3   deg2rad(0);        % outermost tube
         25e-3  deg2rad(0)];     % innermost tube

    % calculate forward kinematics
    arcs = joint2arcparams(tubes, q);   % calcs arc parameters of the deformation
    new_k = 1./(arcs(2:4, 1, 2) * 1e-3);
    new_ang = rad2deg(arcs(3, 2, 2));
    resultant_k(i,:) = new_k';
    resultant_ang(i) = new_ang;
    
%     fprintf("Resultant Radius of Curvature: %.2fmm \n", new_k);
    for j = 1:nTubes
        tubes(j).fwkine(arcs(:,:,j));
    end
end
disp(resultant_k)
disp(resultant_ang)



%% PLOT
plotTubes(tubes);           % regular plot of deformed tubes
plotAllTubes(tubes, q);     % uncomment for plot of each tube undeformed
