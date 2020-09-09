function q = actuator2arcparams(tubes, acts)
%ACTUATOR2ARCPARAMS converts actuator variables of tubes to arc parameters
%  INPUTS
%   tubees = [N]   array of Precurved tube objects
%   acts   = [Nx2] matrix of actuator variables [p theta]
%       p: (m) translation
%       theta: (rad) rotation
%  OUTPUT
%   q  = [2N-1x3xN] matrix of arc parameters for each
%       rows: 2N-1, number of deformed sections (including straight)
%       cols: 3, actuator variables [kappa, theta, trans]
%       sheets: N, each tube

numTubes = length(tubes);
numOverlaps = 2*numTubes-1;

% break down
baseTrans = acts(:,1);
baseRot = acts(:,2);

% get curved section lengths
d = zeros(1, numTubes);
for tube = 1:numTubes
    d(tube) = tubes(tube).Lc;
end

% matrix defining curved/straight/skip for each tube at each overlap
% TODO: make algorithmic
isCurved = zeros(numOverlaps, numTubes);

%% Calculate overlap lengths
% TODO: make algorithmic
% TODO: replace 'isCurved' with function that uses the tube lengths

% init overlap lengths
oLs = zeros(1, numOverlaps);
oLs(1) = baseTrans(2) - baseTrans(1);
oLs(2) = d(1) - oLs(1);

if numTubes == 2
    oLs(3) = d(2) - oLs(2);
    isCurved = [1 0; 1 1; -1 1];
    
elseif numTubes == 3
    oLs(3) = baseTrans(3) - (oLs(2) + oLs(1));
    oLs(4) = d(2) - (oLs(3) + oLs(2));
    oLs(5) = d(3) - oLs(4);
    
    isCurved = [1 0 0; 1 1 0; -1 1 0; -1 1 1; -1 -1 1];
end
    
%% get emergent curvatures for each overlap section
%    each overlap

baseRot
newK = zeros(numOverlaps, 1);  % init new ks for each overlapped section
newTheta = zeros(numOverlaps, 1);
for  link = 1:numOverlaps
    [newK(link) newTheta(link)] = inplane_bending(tubes, baseRot, isCurved(link,:));
end

%% Create q of arc parameters
q = zeros(numOverlaps + 1, 3, numTubes);

for t = 1:numTubes
    
    q(1,:, t) = [0 baseRot(t) baseTrans(1)];    % initial straight section
    
    for link = 1:numOverlaps
        k = newK(link);
        theta = newTheta(link);
        s = oLs(link);
        
        % if the section ends
        if isCurved(link,t) == -1
            k = 0;
            theta = 0;
            s = 0;
        end
            
        q(link+1,:,t) = [k theta s];
    end
end
end

