function q = actuator2arcparams(tubes, acts)
%ACTUATOR2ARCPARAMS converts actuator variables of tubes to arc parameters
%  INPUTS
%   tubees = [N]   array of Precurved tube objects
%   acts   = [Nx2] matrix of actuator variables [p theta]
%       p: (m) translation
%       theta: (rad) rotation
%  OUTPUT
%   q  = [2Nx3xN] matrix of arc parameters for each
%       rows: 2N, number of deformed sections (including initial straight)
%       cols: 3, actuator variables [kappa, theta, trans]
%       sheets: N, each tube

numTubes = length(tubes);
numOverlaps = 2*numTubes;

% break down
baseTrans = acts(:,1);
baseRot = acts(:,2);

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
p = baseTrans; 

% init overlap lengths
Ls = zeros(1, numOverlaps);
Ls(1) = p(1);
Ls(2) = p(2) - p(1);
Ls(3) = d(1) + p(1) - p(2);

if numTubes == 2
    Ls(4) = d(2) - d(1) + p(2) - p(1);
    isCurved = [0 0; 1 0; 1 1; -1 1];
    
elseif numTubes == 3
    Ls(4) = p(3) - p(1) - d(1);
    Ls(5) = d(2) + p(2) - p(3);
    Ls(6) = p(3) + d(3) - p(2) - d(2);
    
    isCurved = [0 0 0; 1 0 0; 1 1 0; -1 1 0; -1 1 1; -1 -1 1];
end
    
%% get emergent curvatures for each overlap section

newK = zeros(numOverlaps, 1);  % init new ks for each overlapped section
newTheta = zeros(numOverlaps, 1);
for  link = 1:numOverlaps
    [newK(link), newTheta(link)] = inplane_bending(tubes, baseRot, isCurved(link,:));
end

newTheta(end+1) = newTheta(end);
newTheta

%% Create q of arc parameters
q = zeros(numOverlaps, 3, numTubes);
for t = 1:numTubes
    
    curr_theta = 0;     % keeps track of current absolute theta
    for link = 1:numOverlaps
        k = newK(link);
        abs_theta = newTheta(link);
        s = Ls(link);
        
        % update relative theta
        if curr_theta ~= abs_theta
            rel_theta = abs_theta - curr_theta;
            curr_theta = curr_theta + rel_theta;
        else
            rel_theta = 0;
        end
        
        % if the section ends
        if isCurved(link,t) == -1 || s == 0
%             disp('skip')
            k = 0;
            s = 0;
        end
        
%         [k rel_theta s]
        q(link,:,t) = [k rel_theta s];
    end
end
end

