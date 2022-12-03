%% COMPARESIMS compares the results of different simulations
%   looks at visibility between simulations in the simulations folder
%   NOTE: must use the same stl model
%
% Author:   Jesse F. d'Almeida <jfdalmeida@wpi.edu>
%         
% Last revision: 7/1/2020

close all, clear, clc

%% OPTIONS
compare = false;

%% Get Simulations
% Ideally would like to compare just 2 simulations, might adjust in future.
% Could also comment this out and hardcode simIDs
simFile = dir('simAnalyzer/simulations/*.mat');
simIDs = convertCharsToStrings({simFile.name});


% Check that there are only 2 simulations
if length(simIDs) ~=2 & compare
    disp("Too many simulations, cannot compare!!!")
    return
end

%% Load in data from simulations

simVars(1) = load(simIDs(1), 'meMesh', 'visibleMapTotal');
for i=2:length(simIDs)
    simVars(i) = load(simIDs(i), 'meMesh', 'visibleMapTotal');
end

% Check if visible map size and stl model
if compare
    for i=1:length(simIDs)-1
        if length(simVars(i).visibleMapTotal) ~= length(simVars(i +1).visibleMapTotal)
            disp("Visible Maps are of different sizes!!!")
            return
        end
        if ~isequal(simVars(i).meMesh, simVars(i+1).meMesh)
            disp("Models are different!!!")
            return
        end
    end
end
meMesh = simVars(1).meMesh;

%% --------------- TOGGLE FOR COMPARISION OF EXACTLY 2 FIGURES --------------
if compare
    % get both visibilities
    map1 = logical(simVars(1).visibleMapTotal);
    map2 = logical(simVars(2).visibleMapTotal)*2;

    % map of faces unique to that map
    unique1 = map1 & ~map2;
    unique2 = (map2 & ~map1)*2;

    % combined of all points
    comp = map1 + map2;

    % create colormap based on required amount
    colors = distinguishable_colors(max(comp)+1);

    % first simulation
    figure('Name', 'Comparing Simulations')
    subplot(131)
    h = stlPlot(meMesh.Vertices * 1e3, meMesh.Faces, simIDs(1), map1, colors);
    set(h, 'DefaultFigureRenderer', 'painters');

    % comparision, can change map to whichever
    subplot(132)
    stlPlot(meMesh.Vertices * 1e3, meMesh.Faces, 'Comparision', unique2, colors);

    % second simulation      
    subplot(133)
    h = stlPlot(meMesh.Vertices * 1e3, meMesh.Faces, simIDs(2), map2, colors);
end
%% ------------ UNCOMMENT FOR MULTIPLE FIGURES DISPLAYED TO BE SAVED AS SVG

colors = [0.2422    0.1504    0.6603;   % purple
          0.9769    0.9839    0.0805];   % yellow
      
for i = 1:length(simIDs)
    map1 = logical(simVars(i).visibleMapTotal);
    meMesh = simVars(i).meMesh;
    
    % plot simulation
    figure('Name', simIDs(i))
    h = stlPlot(meMesh.Vertices * 1e3, meMesh.Faces, simIDs(i), map1, colors);
    set(gca,'FontSize',12);
    set(h, 'DefaultFigureRenderer', 'painters');
    
    %Legend
    unqc = [0 1]';
    label = ["Visible Faces" "Not Visible Faces"];
    hold on;
    for il = 1:length(unqc)
        hl(il) = patch(NaN, NaN, unqc(il));
    end
    
    l = legend(hl,label);
    set(l, 'FontSize', 11);
    set(h, 'DefaultFigureRenderer', 'painters');
end



