
function [robot, k, kj, theta, thetaj, s] = get_Curvature(p)
    % stainless steel 304 material properties
    yield = 290e6;
    modulus = 193e9;
    
    
    % pla
%     yield = 43e6;
%     modulus = 4.1e9;
%     emax = yield/modulus;   % strain at yield
    % Aluminum
%     yield = 276e6;
%     modulus = 68.9e9;
    
%     emax = yield/modulus;   % strain at yield
        emax = p
    % nitinol
%     yield = 200e6;
%     modulus = 40e9;
%     emax = 0.08;

    % tube params
    OD = 2e-3;
    ID = 1.6e-3;
    ro = OD/2;      % radius of outer wall
    ri = ID/2;      % radius of inner wall
    n = 15;          % number of cutouts
    g = OD * .79;      % depth of cut based on percent of diameter cut
    u = 5 * 1e-3;   % height of uncut section
    h = 5 * 1e-3;   % height of cut section
    
    % assemble wrist
    cutouts.w = g * ones(1,n); % [m]
    cutouts.u = u * ones(1,n); % [m]
    cutouts.h = h * ones(1,n); % [m]
    cutouts.alpha = zeros(1,n);
    robot = Wrist(OD, ID, n, cutouts);
    
    ybar = robot.ybar(1);   % all ybars are the same for each notch
    
    % location of max strain on cut section
    % max dist away from ybar is either the inner or outer wall 
    y = [ro - ybar, ybar - (g - ro)];
    [d, idx] = max(abs(y));   
    d = d * sign(y(idx));
    
    kj = emax/(d - ybar * emax);    % calc max kappa
    
    thetaj = (h/(1 + ybar * kj))*kj;    % bending angle of each notch
    theta = thetaj * n;                 % total bending angle of tube
    
    sj = h / (1 + ybar * kj);    % arc length of each notch
    s = sj * n + u * (n-1);
    
    % calc min row
    rho = s/theta;
    k = 1/rho;
    
end
