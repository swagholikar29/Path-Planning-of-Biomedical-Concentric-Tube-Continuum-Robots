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
p = acts(:,1);
theta = acts(:,2);

% get curved section lengths
d = zeros(1, numTubes);
for i = 1:numTubes
    d(i) = tubes(i).Lc;
end

% matrix defining curved/straight/skip for each tube at each overlap
% TODO: make algorithmic
isCurved = zeros(numOverlaps, numTubes);

%% Calculate overlap lengths
% TODO: make algorithmic

% init overlap lengths
oLs = zeros(1, numOverlaps);
oLs(1) = p(2) - p(1);
oLs(2) = d(1) - oLs(1);

if numTubes == 2
    oLs(3) = d(2) - oLs(2);
    isCurved = [1 0; 1 1; -1 1];
    
elseif numTubes == 3
    oLs(3) = p(3) - (oLs(2) + oLs(1));
    oLs(4) = d(2) - (oLs(3) + oLs(2));
    oLs(5) = d(3) - oLs(4);
    
    isCurved = [1 0 0; 1 1 0; -1 1 0; -1 1 1; -1 -1 1];
end
    
%% get emergent curvatures for each overlap section
%    each overlap

newK = zeros(numOverlaps, 1);  % init new ks for each overlapped section
for  i = 1:numOverlaps
    newK(i) = inplane_bending(tubes, isCurved(i,:));
end

%% Create q of arc parameters
q = zeros(numOverlaps + 1, 3, numTubes);

for i = 1:numTubes
    q(:,:,i) = deformation(newK, oLs, p(1), isCurved(:,i)); 
end
end

