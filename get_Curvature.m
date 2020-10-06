
function [robot, k, kj, theta, thetaj, s] = get_Curvature(p)
    % stainless steel 304 material properties
    yield = 290e6;
    modulus = 193e9;
    
    % nitinol
%     yield = 200e6;
%     modulus = 40e9;
    
    % pla
%     yield = 43e6;
%     modulus = 4.1e9;
    emax = yield/modulus;   % strain at yield
    
    % tube params
    OD = 3e-3;
    ID = 2.8e-3;
    ro = OD/2;  % radius of outer wall
    ri = ID/2;  % radius of inner wall
    n = 20;      % number of cutouts
    g = OD*p;   % depth of cut based on percent of diameter cut
    u = 2e-3;   % height of uncut section
    h = 2e-3;   % height of cut section
    
    % assemble wrist
    cutouts.w = g * ones(1,n); % [m]
    cutouts.u = u * ones(1,n); % [m]
    cutouts.h = h * ones(1,n); % [m]
    cutouts.alpha = zeros(1,n);
    robot = Wrist(OD, ID, n, cutouts);
    
    ybar = robot.ybar(1);   % all ybars are the same for each notch
    
    
    d = max([g - ro - ybar, ro - ybar]);    % location of max strain on cut section
    
    kj = emax/(d - ybar* emax);   % calc max kappa
    
    thetaj = (h/(1 + ybar * kj))*kj;      % bending angle of each notch
    theta = thetaj * n;     %  total bending angle of tube
    s = h*n + u*(n-1);      % total arc length of tube 
   
    k = theta/s;
    
end
