function arcs = joint2arcparams(tubes, q)
%ACTUATOR2ARCPARAMS converts actuator variables of tubes to arc parameters
%  INPUTS
%   tubees = [N]   array of Precurved tube objects
%   q      = [Nx2] matrix of actuator variables [p theta]
%       p: (m) translation
%       alpha: (rad) rotation
%  OUTPUT
%   q  = [2Nx3xN] matrix of arc parameters for each
%       rows: 2N, number of deformed sections (including initial straight)
%       cols: 3, actuator variables [kappa, phi, s]
%       sheets: N, each tube

numTubes = length(tubes);
numOverlaps = 2*numTubes;

% break down
p = q(:,1);
alpha = q(:,2);

% get curved section lengths
d = zeros(1, numTubes);
for tube = 1:numTubes
    d(tube) = tubes(tube).Lc;
end

% matrix defining curved/straight/skip for each tube at each overlap
isCurved = zeros(numOverlaps, numTubes);

%% Calculate overlap lengths
% TODO: make algorithmic
% TODO: replace 'isCurved' with function that uses the tube lengths to
%       calculate

[Ls isCurved] = calcLinkLengths(tubes, p)

% init overlap lengths
% Ls = zeros(1, numOverlaps);
% Ls(1) = p(1);
% Ls(2) = p(2) - p(1);
% Ls(3) = d(1) + p(1) - p(2);
% 
% if numTubes == 2
%     Ls(4) = d(2) - d(1) + p(2) - p(1);
%     isCurved = [0 0; 1 0; 1 1; -1 1];
%     
% elseif numTubes == 3
%     Ls(4) = p(3) - p(1) - d(1);
%     Ls(5) = d(2) + p(2) - p(3);
%     Ls(6) = p(3) + d(3) - p(2) - d(2);
%     
%     isCurved = [0 0 0; 1 0 0; 1 1 0; -1 1 0; -1 1 1; -1 -1 1];
% end
    
%% get emergent curvatures for each overlap section

newK = zeros(numOverlaps, 1);  % init new ks for each overlapped section
phi = zeros(numOverlaps, 1);
for  link = 1:numOverlaps
    [newK(link), phi(link)] = inplane_bending(tubes, alpha, isCurved(link,:));
end

%% Create arc parameters
arcs = zeros(numOverlaps, 3, numTubes);
for t = 1:numTubes
    
    curr_phi = 0;     % keeps track of current angle
    for link = 1:numOverlaps
        k = newK(link);
        abs_phi = phi(link);
        s = Ls(link);
        
        % update relative theta
        if curr_phi ~= abs_phi
            rel_phi = abs_phi - curr_phi;
            curr_phi = curr_phi + rel_phi;
        else
            rel_phi = 0;
        end
        
        % if the section ends
        if isCurved(link,t) == -1 || s == 0
            k = 0;
            s = 0;
        end
        
        arcs(link,:,t) = [k rel_phi s];
    end
end
end

