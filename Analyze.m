%% Script to create series of precurved concentric tubes and evaluate their kinematics and mechanics
close all, clear, clc
addpath('kinematics');

%% CREATE PRECURVED TUBES
% get resultant curvatures from list of tubes
Rs = [68.53	44.53];


Ks = 1 ./ (Rs * 1e-3);         % (1/m) curvatures
straight = zeros(length(Ks),1);

% Ks = [straight Ks];

ODs = [3.8e-3 2.2e-3];         % (m) outer diameters 

numRs = size(Ks,1);
resultant_k = zeros(size(Ks,1),1);

for i = 1:numRs
    % parameters for generating tubes
    nTubes = 2;                     % number of tubes
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
         0e-3  deg2rad(0)];     % innermost tube

    % calculate forward kinematics
    arcs = joint2arcparams(tubes, q);   % calcs arc parameters of the deformation
    new_k = 1/(arcs(3, 1, 1) *1e-3);
    resultant_k(i) = new_k;
    
%     fprintf("Resultant Radius of Curvature: %.2fmm \n", new_k);
    for j = 1:nTubes
        tubes(j).fwkine(arcs(:,:,j));
    end
end
disp(resultant_k)



%% PLOT
% plotTubes(tubes);           % regular plot of deformed tubes
% plotAllTubes(tubes, q);     % uncomment for plot of each tube undeformed
