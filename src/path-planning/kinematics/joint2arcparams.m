function arcs = joint2arcparams(tubes, q)
%ACTUATOR2ARCPARAMS converts actuator variables of tubes to arc parameters
%  INPUTS
%   tubes = [N]   array of Precurved tube objects
%   q     = [Nx2] matrix of actuator variables [p theta]
%   p     : (m) translation
%   alpha : (rad) rotation
%  OUTPUT
%   q  = [2Nx3xN] matrix of arc parameters for each
%       rows: 2N, number of deformed sections (including initial straight)
%       cols: 3, actuator variables [kappa, phi, s]
%       sheets: N, each tube

numTubes = length(tubes);
numOverlaps = 2*numTubes; % N= 2n 

% break down
p = q(:,1);
alpha = q(:,2);

% get curved section lengths
d = zeros(1, numTubes);
for tube = 1:numTubes
    d(tube) = tubes(tube).Lc; %Lc = length of curved section
end

% matrix defining curved/straight/skip for each tube at each overlap
isCurved = zeros(numOverlaps, numTubes);

%% Calculate overlap lengths

[Ls isCurved] = calcLinkLengths(tubes, p);
    
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

